## Functions for glycopeptide analysis using MSFragger
## Bradberry et al., 2023: "N-glycoproteomics of brain synapses and synaptic vesicles"
## These scripts are intended for use with functions contained in an accompanying R script file
## Please post to Github or e-mail mazdak dot bradberry at nyspi.columbia.edu with any questions


# Load label-free quantification data and generate intensity lists (Fig. 1)
LFQ_SV <- read.delim("SV_summary.csv", sep=",", stringsAsFactors = FALSE)
LFQ_syn <- read.delim("synapt_summary.csv", sep=",", stringsAsFactors = FALSE)
SVintensityList <- data.frame(Entry.Name = LFQ_SV$Entry, Intensity = LFQ_SV$Intensity.average)
Synintensitylist <- data.frame(Entry.Name = LFQ_syn$Entry, Intensity = LFQ_syn$Intensity.average)

# Import a table with organelle localization protein entries (Fig. 2)
organelle_entries <- read.delim("Sigma_organelles.csv", sep=",", header = TRUE)

# Create a table with log(10) LFQ values for SVs and synaptosomes (Fig. 2)
x <- merge(LFQ_SV, LFQ_syn, by = "Entry.Name", all = TRUE)
LFQ_for_organelles <- data.frame(Entry.Name = x$Entry.Name, 
                                 SV_1 = log10(x$SV_1.Intensity),
                                 SV_2 = log10(x$SV_2.Intensity),
                                 SV_3 = log10(x$SV_3.Intensity),
                                 syn_1 = log10(x$synapt_1.Intensity),
                                 syn_2 = log10(x$synapt_2.Intensity),
                                 syn_3 = log10(x$synapt_3.Intensity)
)

# Add LFQ values to the organelle localization table (Fig. 2)
all_organelles <- merge(organelle_entries, LFQ_for_organelles, by = "Entry.Name", all.x = TRUE)
write.csv(all_organelles, "all_organelles.csv")


# Compare LFQ data between the SV2-IP samples in this study and Syt1-IP samples from Bradberry 2022 (Fig. S1)
Syt1oldintensitylist <- read.delim("Bradberry_Table_5-1.csv", sep = ",")
Syt1oldintensitylist <- subset(Syt1oldintensitylist, Syt1_1.Intensity*Syt1_2.Intensity > 0)
SV2_nonzero <- subset(LFQ_SV, SV_1.Intensity*SV_2.Intensity*SV_3.Intensity > 0)
SVIP_compare <- merge(SV2_nonzero, Syt1oldintensitylist, by = "Entry.Name")
SVIP_compare$SV2_log10 <- log10(SVIP_compare$Intensity.average)
SVIP_compare$Syt1_log10 <- log10(SVIP_compare$Syt1_average)
write.csv(SVIP_compare, "SVIP_compare.csv")

Syt1_detected <- subset(Syt1oldintensitylist, Syt1_average > 0)
SV_detected <- subset(LFQ_SV, Intensity.average > 0)
nrow(Syt1_detected)
nrow(SV_detected)
sum(Syt1_detected$Entry.Name %in% SV_detected$Entry.Name)


# Load glycoPSM tables (Found in Supplemental Table S2)
SV_E<-read.delim("psm_SV_enriched.tsv", sep="\t", stringsAsFactors = FALSE)
SV_NE<-read.delim("psm_SV_NE.tsv", sep="\t", stringsAsFactors = FALSE)
Syn_sce<-read.delim("psm_sce_120_1.tsv", sep="\t", stringsAsFactors = FALSE)
Syn_AIETD<-read.delim("psm_AIETD_120_2.tsv", sep="\t", stringsAsFactors = FALSE)
syt1_SV2<-read.delim("psm_syt1_SV2.tsv", sep = "\t", stringsAsFactors = FALSE)
Syn<-rbind(Syn_AIETD, Syn_sce)


# Load glycan database (format: column 1 = MSfragger glycan ID, column 2 = annotated glycan type) for annotation (Supplemental Table S4)
database <- read.delim("glycanList_annotated_simple.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)


# Process glycoPSMs and generate output (Figs. 4, 6; Supplemental Tables S3, S5, S6)
SV_E_clean <- psmCleanup(SV_E, 0.025, TRUE)
SV_E_unique <- uniqueGlycos(SV_E_clean)
SV_E_included <- glycanAnnotate(SV_E_unique, TRUE, TRUE)
SV_E_included <- intensityAssign(SV_E_included, SVintensityList)
SV_E_ntile <- ntile(SVintensityList, SV_E_included, 5)
SV_E_included$percentile <- sapply(SV_E_included$Intensity, ntileAssign, SV_E_ntile)
SV_E_occurenceTable <- glycoCountPercentile(SV_E_included)

SV_E_all_annotated <- glycanAnnotate(SV_E_clean, FALSE, TRUE)
SV_E_all_table <- data.frame(table(SV_E_all_annotated$Annotated))
SV_E_unique_annotated <- glycanAnnotate(SV_E_unique, FALSE, TRUE)
SV_E_unique_table <-data.frame(table(SV_E_unique_annotated$Annotated))


SV_NE_clean<-psmCleanup(SV_NE, 0.025, TRUE)
SV_NE_unique <- uniqueGlycos(SV_NE_clean)
SV_NE_included <- glycanAnnotate(SV_NE_unique, TRUE, TRUE)
SV_NE_included <- intensityAssign(SV_NE_included, SVintensityList)
SV_NE_ntile <- ntile(SVintensityList, SV_NE_included, 5)
SV_NE_included$percentile <- sapply(SV_NE_included$Intensity, ntileAssign, SV_NE_ntile)
SV_NE_occurenceTable <- glycoCountPercentile(SV_NE_included)

SV_NE_all_annotated <- glycanAnnotate(SV_NE_clean, FALSE, TRUE)
SV_NE_all_table <- data.frame(table(SV_NE_all_annotated$Annotated))
SV_NE_unique_annotated <- glycanAnnotate(SV_NE_unique, FALSE, TRUE)
SV_NE_unique_table <-data.frame(table(SV_NE_unique_annotated$Annotated))


Syn<-rbind(Syn_AIETD, Syn_sce)
Syn_clean<-psmCleanup(Syn, 0.025, TRUE)
Syn_unique <- uniqueGlycos(Syn_clean)
Syn_included <- glycanAnnotate(Syn_unique, TRUE, TRUE)
Syn_included <- intensityAssign(Syn_included, Synintensitylist)
Syn_ntile <- ntile(Synintensitylist, Syn_included, 5)
Syn_included$percentile <- sapply(Syn_included$Intensity, ntileAssign, Syn_ntile)
Syn_occurenceTable <- glycoCountPercentile(Syn_included)

Syn_all_annotated <- glycanAnnotate(Syn_clean, FALSE, TRUE)
Syn_all_table <- data.frame(table(Syn_all_annotated$Annotated))
Syn_unique_annotated <- glycanAnnotate(Syn_unique, FALSE, TRUE)
Syn_unique_table <-data.frame(table(Syn_unique_annotated$Annotated))

syt1_SV2_clean <- psmCleanup(syt1_SV2, 0.025, TRUE)
syt1_SV2_annotated <- glycanAnnotate(syt1_SV2_clean, FALSE, TRUE)

all_psms <- rbind(SV_E_clean, SV_NE_clean, Syn_clean, syt1_SV2_clean)
all_unique <- uniqueGlycos(all_psms)
all_annotated <- glycanAnnotate(all_unique, FALSE, TRUE)
selected_Unique <- subset(all_annotated, Entry.Name %in% c("SV2A_MOUSE", "SV2B_MOUSE", 
                                                           "SV2C_MOUSE","SYPH_MOUSE",
                                                           "SYT1_MOUSE", "THY1_MOUSE",
                                                           "SYNPR_MOUSE"))

write.csv(apply(SV_E_included, 2, as.character), "SV_E_included.csv")
write.csv(apply(SV_NE_included, 2, as.character), "SV_NE_included.csv")
write.csv(apply(Syn_included, 2, as.character), "Syn_included.csv")

write.csv(apply(SV_E_all_annotated, 2, as.character), "SV_E_all_annotated.csv")
write.csv(apply(SV_NE_all_annotated, 2, as.character), "SV_NE_all_annotated.csv")
write.csv(apply(Syn_all_annotated, 2, as.character), "Syn_all_annotated.csv")
write.csv(apply(syt1_SV2_annotated, 2, as.character), "Syt1_SV_IP_all_annotated.csv")

write.csv(SV_E_occurenceTable, "SV_E_occurenceTable.csv")
write.csv(SV_NE_occurenceTable, "SV_NE_occurenceTable.csv")
write.csv(Syn_occurenceTable, "Syn_occurenceTable.csv")

write.csv(SV_E_unique_table, "SV_E_unique_table.csv")
write.csv(SV_NE_unique_table, "SV_NE_unique_table.csv")
write.csv(Syn_unique_table, "Syn_unique_table.csv")

write.csv(SV_E_all_table, "SV_E_all_table.csv")
write.csv(SV_NE_all_table, "SV_NE_all_table.csv")
write.csv(Syn_all_table, "Syn_all_table.csv")

write.csv(selected_Unique, "Selected_unique.csv")


# Generate tables of all unique glycoPSMs (Supplementary Figs. S2, S3)
SV <- rbind(SV_E, SV_NE)
SV <- psmCleanup(SV, 0.025, TRUE)
SV <- glycanAnnotate(SV, FALSE, TRUE)
SV_uniqueGlycos <- subset(SV, !duplicated(Glycan))
SV_table <- data.frame(table(SV$Annotated))
write.csv(SV_uniqueGlycos, "SV_labeled.csv")
write.csv(SV_table, "SV_table.csv")

Syn <-rbind(Syn_AIETD, Syn_sce)
Syn <- psmCleanup(Syn, 0.025, TRUE)
Syn <- glycanAnnotate(Syn, FALSE, TRUE)
Syn_uniqueGlycos <- subset(Syn, !duplicated(Glycan))
Syn_table<-data.frame(table(Syn$Annotated))
write.csv(Syn_uniqueGlycos, "Syn_labeled.csv")
write.csv(Syn_table, "Syn_table.csv")


# Define the number of unique glycoproteins and glycopeptides
length(unique(c(SV$Entry.Name, SV_NE$Entry.Name, Syn$Entry.Name)))
length(unique(c(SV_unique$Glycopeptide, SV_NE_unique$Glycopeptide, Syn_unique$Glycopeptide)))


# Compare overlap of glycopeptides in synaptosome and SV samples (Fig. 7)
SV <- rbind(SV_E, SV_NE)
SV <- psmCleanup(SV, 0.025, TRUE)
SV_unique <- uniqueGlycos(SV)
sum(Syn_unique$Glycopeptide %in% SV_unique$Glycopeptide)


# Perform glycopeptide and glycoprotein overlap analysis by mannosylation/fucosylation (Fig. 7)
all_annotated_hifuc <- subset(all_annotated, Annotated %in% c("Fuc2", "Fuc3"))
all_annotated_onefuc <- subset(all_annotated, Annotated %in% c("Fuc1"))
all_annotated_man <-subset (all_annotated, Annotated %in% c("Man5", "Man6", "Man7", "Man8", "Man9"))

glycoproteins_hifuc <- unique(all_annotated_hifuc$Entry.Name)
glycoproteins_onefuc <- unique(all_annotated_onefuc$Entry.Name)
glycoproteins_man <- unique(all_annotated_man$Entry.Name)
nonfuc <- subset(glycoproteins_man, ! glycoproteins_man %in% c(glycoproteins_hifuc, glycoproteins_onefuc))


onefuc_hifuc <- subset(glycoproteins_onefuc, glycoproteins_onefuc %in% glycoproteins_hifuc)
onefuc_man <-subset(glycoproteins_onefuc, glycoproteins_onefuc %in% glycoproteins_man)
hifuc_man <- subset(glycoproteins_hifuc, glycoproteins_hifuc %in% glycoproteins_man)

sum(onefuc_hifuc %in% glycoproteins_man)
sum(glycoproteins_onefuc %in% glycoproteins_hifuc)

sum(glycoproteins_onefuc %in% glycoproteins_man)
sum(glycoproteins_hifuc %in% glycoproteins_man)


# Write glycoprotein lists for GO analysis (Fig. 7)
write.csv(glycoproteins_hifuc, "glycoproteins_hifuc.csv")
write.csv(glycoproteins_onefuc, "glycoproteins_onefuc.csv")
write.csv(glycoproteins_man, "glycoproteins_man.csv")
write.csv(nonfuc, "glycoproteins_nonfuc.csv")


