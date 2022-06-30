## Load glycan database for annotation
database <- read.delim("glycanList_annotated_noNeu.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
## Remove unneeded columns from PSM table and filter PSMs by glycan Q-value
## Add a "Glycopeptide" column that contains both the peptide sequence and glycan composition
psmCleanup <- function(psm, q, removeAbs) {
  cleanpsm <- data.frame(Peptide <- psm$Peptide)
  colnames(cleanpsm) <- "Peptide"
  cleanpsm$Glycan <- psm$Observed.Modifications
  cleanpsm$Glycan.Q <- psm$Glycan.q.value
  cleanpsm$Glycan.score <- psm$Glycan.Score
  cleanpsm$Entry.Name <- psm$Entry.Name
  cleanpsm$Glycopeptide <- paste(psm$Peptide, cleanpsm$Glycan, sep="_")
  cleanpsm <- subset(cleanpsm, Glycan.Q < q)
  if (removeAbs == TRUE) {
    cleanpsm <- subset(cleanpsm, ! Entry.Name %in% c("IGHG1_MOUSE",
                                                     "IGG2B_MOUSE",
                                                     "IGHM_MOUSE",
                                                     "IGHG3_MOUSE"))
  }
  return(cleanpsm)
}
## Removes duplicated glycopeptides from a cleaned-up PSM table
uniqueGlycos <- function(psm) {
  unique <- subset(psm, !duplicated(Glycopeptide))
  return(unique)
}
## Function to apply to glycan composition to assign annotation
glycanID <- function(comp, database) {
  return(as.character(database$id[comp == database$comp]))
}
# Applies glycanID function to psm table, generating a column with glycan annotation
# Removes non-annotated glycans
# Option to remove antibody-derived glycopeptides
glycanAnnotate <- function(psm, removeAbs) {
  psm$Annotated <- as.character(sapply(psm$Glycan, "glycanID", database = database))
  annotatedGlyco <- subset(psm, Annotated != "character(0)")
  if (removeAbs == TRUE) {
    annotatedGlyco <- subset(annotatedGlyco, ! Entry.Name %in% c("IGHG1_MOUSE",
                                                                 "IGG2B_MOUSE",
                                                                 "IGHM_MOUSE",
                                                                 "IGHG3_MOUSE"))
  }
  return(annotatedGlyco)
}
# Generates intensities corresponding to n intensity percentile groups using a list of protein intensities
ntile <- function(intensityList, annotatedGlyco, n) {
  ntileList <- subset(intensityList, intensityList$Entry.Name %in% annotatedGlyco$Entry.Name)
  GPntile <- quantile(ntileList$Intensity, probs = seq(0, 1, by = 1/n), names = TRUE)
  return (GPntile)
}
# Imports intensity values into annotated glycoPSM table
intensityAssign <- function(intensityList, annotatedGlyco) {
  GPintensity <- subset(intensityList, Entry.Name %in% annotatedGlyco$Entry.Name)
  GPintensity <- merge(annotatedGlyco, GPintensity)
  return(GPintensity)
}
# Function to look up the percentile for a given intensity value
ntileAssign <- function(intensity, GPntile) {
  for (i in 2:length(GPntile)) {
    if ((intensity < GPntile[i]) & (intensity > GPntile[i - 1])) {
      return(names(GPntile[i-1]))
    }
  }
}
# Function to normalize a numeric vector to the sum of its components
normalize <- function (vec) {
  newVec = vec/sum(vec)
  return(newVec)
}
# Generates a table that tabulates all annotated glycans in an annonated PSM table
# And also tabulates the annotated glycans by percentile group
# Then cleans up the output table and generates normalized values for each percentile
glycoCountPercentile <- function(annotatedGlyco) {
  occurenceTable <- data.frame(table(annotatedGlyco$Annotated))
  occurenceTable <- setNames(occurenceTable, c("AnnotatedGlycan", "Frequency"))
  for (i in unique(annotatedGlyco$percentile)) {
    ntileTable <- subset(annotatedGlyco, percentile == i)
    ntileTable <- data.frame(table(ntileTable$Annotated))
    ntileTable <- setNames(ntileTable, c("AnnotatedGlycan", paste(i)))
    occurenceTable <- merge(occurenceTable, ntileTable, all.x = TRUE)
  }
  occurenceTable[is.na(occurenceTable)] = 0
  rownames(occurenceTable) <- occurenceTable$AnnotatedGlycan
  occurenceTable <- subset(occurenceTable, select = -AnnotatedGlycan)
  occurenceTable <- subset(occurenceTable, rownames(occurenceTable) %in% c("Man5", "Man6", "Man7", "Man8", "Man9", "Fuc0", "Fuc1", "Fuc2", "Fuc3"))
  glycoNorm <- data.frame(sapply(occurenceTable, normalize))
  rownames(glycoNorm) <- paste("norm",rownames(occurenceTable))
  colnames(glycoNorm) <- colnames(occurenceTable)
  occurenceTable <- rbind(occurenceTable,glycoNorm)
  return(occurenceTable)
}