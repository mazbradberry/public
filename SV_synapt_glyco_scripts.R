
# Load label-free quantification data
LFQ_SV <- read.delim("SV_summary.csv", sep=",", stringsAsFactors = FALSE)
LFQ_syn <- read.delim("synapt_summary.csv", sep=",", stringsAsFactors = FALSE)
SVintensityList <- data.frame(Entry.Name = LFQ_SV$Entry, Intensity = LFQ_SV$Intensity.average)
Synintensitylist <- data.frame(Entry.Name = LFQ_syn$Entry, Intensity = LFQ_syn$Intensity.average)

# Load glycoPSM tables
SV<-read.delim("psm_SV_enriched.tsv", sep="\t", stringsAsFactors = FALSE)
SV_NE<-read.delim("psm_SV_NE.tsv", sep="\t", stringsAsFactors = FALSE)
Syn_sce<-read.delim("psm_sce_120_1.tsv", sep="\t", stringsAsFactors = FALSE)
Syn_AIETD<-read.delim("psm_AIETD_120_2.tsv", sep="\t", stringsAsFactors = FALSE)
syt1_SV2<-read.delim("psm_syt1_SV2.tsv", sep = "\t", stringsAsFactors = FALSE)
Syn<-rbind(Syn_AIETD, Syn_sce)

# Process glycoPSMs and generate output
SV <- psmCleanup(SV, 0.025, TRUE)
SV_unique <- uniqueGlycos(SV)
SV_annotated <- glycanAnnotate(SV_unique, TRUE)
SV_annotated <- intensityAssign(SV_annotated, SVintensityList)
SV_ntile <- ntile(SVintensityList, SV_annotated, 5)
SV_annotated$percentile <- sapply(SV_annotated$Intensity, ntileAssign, SV_ntile)
SV_occurenceTable <- glycoCountPercentile(SV_annotated)
SV_all_annotated <- glycanAnnotate(SV, TRUE)
SV_all_table <- data.frame(table(SV_all_annotated$Annotated))

SV_NE<-psmCleanup(SV_NE, 0.025, TRUE)
SV_NE_unique <- uniqueGlycos(SV_NE)
SV_NE_annotated <- glycanAnnotate(SV_NE_unique, TRUE)
SV_NE_annotated <- intensityAssign(SV_NE_annotated, SVintensityList)
SV_NE_ntile <- ntile(SVintensityList, SV_NE_annotated, 5)
SV_NE_annotated$percentile <- sapply(SV_NE_annotated$Intensity, ntileAssign, SV_NE_ntile)
SV_NE_occurenceTable <- glycoCountPercentile(SV_NE_annotated)
SV_NE_all_annotated <- glycanAnnotate(SV_NE, TRUE)
SV_NE_all_table <- data.frame(table(SV_NE_all_annotated$Annotated))

Syn<-psmCleanup(Syn, 0.025, TRUE)
Syn_unique <- uniqueGlycos(Syn)
Syn_annotated <- glycanAnnotate(Syn_unique, TRUE)
Syn_annotated <- intensityAssign(Syn_annotated, Synintensitylist)
Syn_ntile <- ntile(Synintensitylist, Syn_annotated, 5)
Syn_annotated$percentile <- sapply(Syn_annotated$Intensity, ntileAssign, Syn_ntile)
Syn_occurenceTable <- glycoCountPercentile(Syn_annotated)
Syn_all_annotated <- glycanAnnotate(Syn, TRUE)
Syn_all_table <- data.frame(table(Syn_all_annotated$Annotated))

syt1_SV2 <- psmCleanup(syt1_SV2, 0.025, TRUE)
all_psms <- rbind(SV, SV_NE, Syn, syt1_SV2)
all_unique <- uniqueGlycos(all_psms)
all_annotated <- glycanAnnotate(all_unique, TRUE)
selected_Unique <- subset(all_annotated, Entry.Name %in% c("SV2A_MOUSE", "SV2B_MOUSE", 
                                                           "SV2C_MOUSE","SYPH_MOUSE",
                                                           "SYT1_MOUSE", "THY1_MOUSE",
                                                           "SYNPR_MOUSE"))


write.csv(SV_occurenceTable, "SV_occurenceTable.csv")
write.csv(SV_NE_occurenceTable, "SV_NE_occurenceTable.csv")
write.csv(Syn_occurenceTable, "Syn_occurenceTable.csv")

write.csv(SV_all_table, "SV_all_table.csv")
write.csv(SV_NE_all_table, "SV_NE_all_table.csv")
write.csv(Syn_all_table, "Syn_all_table.csv")

write.csv(selected_Unique, "Selected_unique.csv")


# Define the number of unique glycoproteins and glycopeptides
length(unique(c(SV$Entry.Name, SV_NE$Entry.Name, Syn$Entry.Name)))
length(unique(c(SV_unique$Glycopeptide, SV_NE_unique$Glycopeptide, Syn_unique$Glycopeptide)))



# Perform glycopeptide and glycoprotein overlap analysis

SV_pooled <- rbind(SV, SV_NE)
SV_pooled_unique <- uniqueGlycos(SV_pooled)
sum(Syn_unique$Glycopeptide %in% SV_pooled_unique$Glycopeptide)

all_annotated_hifuc <- subset(all_annotated, Annotated %in% c("Fuc2", "Fuc3"))
all_annotated_onefuc <- subset(all_annotated, Annotated %in% c("Fuc1"))
all_annotated_man <-subset (all_annotated, Annotated %in% c("Man5", "Man6", "Man7", "Man8", "Man9"))

glycoproteins_hifuc <- unique(all_annotated_hifuc$Entry.Name)
glycoproteins_onefuc <- unique(all_annotated_onefuc$Entry.Name)
glycoproteins_man <- unique(all_annotated_man$Entry.Name)
nonfuc <- subset(glycoproteins_man, ! glycoproteins_man %in% onefuc_hifuc)

onefuc_hifuc <- subset(glycoproteins_onefuc, glycoproteins_onefuc %in% glycoproteins_hifuc)
onefuc_man <-subset(glycoproteins_onefuc, glycoproteins_onefuc %in% glycoproteins_man)
hifuc_man <- subset(glycoproteins_hifuc, glycoproteins_hifuc %in% glycoproteins_man)

sum(onefuc_hifuc %in% glycoproteins_man)
sum(glycoproteins_onefuc %in% glycoproteins_hifuc)

sum(glycoproteins_onefuc %in% glycoproteins_man)
sum(glycoproteins_hifuc %in% glycoproteins_man)


# Write glycoprotein lists for GO analysis

write.csv(glycoproteins_hifuc, "glycoproteins_hifuc.csv")
write.csv(glycoproteins_onefuc, "glycoproteins_onefuc.csv")
write.csv(glycoproteins_man, "glycoproteins_man.csv")
write.csv(nonfuc, "glycoproteins_nonfuc.csv")


# Filter LFQ results by glycoproteins depending on fucosylation

LFQ_syn_hifuc <- subset(LFQ_syn, Entry.Name %in% glycoproteins_hifuc)
LFQ_SV_hifuc <- subset(LFQ_SV, Entry.Name %in% glycoproteins_hifuc)

write.csv(LFQ_syn_hifuc, "LFQ_syn_hifuc.csv")
write.csv(LFQ_SV_hifuc, "LFQ_SV_hifuc.csv")




