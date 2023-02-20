## Functions for glycopeptide analysis using MSFragger
## Bradberry et al., 2023: "N-glycoproteomics of brain synapses and synaptic vesicles"
## These functions are intended for use with scripts contained in an accompanying R script file
## Please post to Github or e-mail mazdak dot bradberry at nyspi.columbia.edu with any questions


# Generate a cleaned-up PSM table
# Remove unneeded columns from MSFragger psm table and filter PSMs by glycan Q-value
# Add a "Glycopeptide" column that contains the concatenated peptide sequence and glycan composition
# Specify in input: q-value cutoff (q) and whether to remove antibody-derived PSMs (removeAbs)
psmCleanup <- function(psm, q, removeAbs) {
  cleanpsm <- data.frame(Peptide <- psm$Peptide)
  colnames(cleanpsm) <- "Peptide"
  cleanpsm$Glycan <- psm$Observed.Modifications
  cleanpsm$Glycan.Q <- psm$Glycan.q.value
  cleanpsm$Glycan.score <- psm$Glycan.Score
  cleanpsm$Entry.Name <- psm$Entry.Name
  cleanpsm$Entry.Description <- psm$Protein.Description
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
# Remove duplicated glycopeptides from a cleaned-up PSM table
uniqueGlycos <- function(psm) {
  unique <- subset(psm, !duplicated(Glycopeptide))
  return(unique)
}

# Assign annotation to each glycan composition based on the glycan database
glycanID <- function(comp, database) {
  return(as.character(database$id[comp == database$comp]))
}

# Apply glycanID function to a cleaned-up psm table, generating a column with glycan annotation
# Specify in input whether to remove non-annotated glycoPSMs from the table (removeNonAnn)
# Specify in input whether to remove antibody-derived PSMs (removeAbs)
glycanAnnotate <- function(psm, removeNonAnn, removeAbs) {
  psm$Annotated <- as.character(sapply(psm$Glycan, "glycanID", database = database))
  annotatedGlyco <- psm
  if (removeNonAnn == TRUE) {
    annotatedGlyco <- subset(psm, Annotated != "character(0)")
  }
  if (removeAbs == TRUE) {
    annotatedGlyco <- subset(annotatedGlyco, ! Entry.Name %in% c("IGHG1_MOUSE",
                                                                 "IGG2B_MOUSE",
                                                                 "IGHM_MOUSE",
                                                                 "IGHG3_MOUSE"))
  }
  return(annotatedGlyco)
}


# Generate intensity quantile groups for a set of glycoproteins from a list of protein intensities 
# Specify in input table containing glycoproteins (annotatedGlyco) and how many quantile groups to include (n)
ntile <- function(intensityList, annotatedGlyco, n) {
  ntileList <- subset(intensityList, intensityList$Entry.Name %in% annotatedGlyco$Entry.Name)
  GPntile <- quantile(ntileList$Intensity, probs = seq(0, 1, by = 1/n), names = TRUE)
  return (GPntile)
}


# Import intensity values from a protein intensity list into annotated glycoPSM table
intensityAssign <- function(intensityList, annotatedGlyco) {
  GPintensity <- subset(intensityList, Entry.Name %in% annotatedGlyco$Entry.Name)
  GPintensity <- merge(annotatedGlyco, GPintensity)
  return(GPintensity)
}


# Look up the intensity quantile for a given glycoprotein intensity value
ntileAssign <- function(intensity, GPntile) {
  for (i in 1:(length(GPntile)-1) ) {
    if ((intensity <= GPntile[i+1]) & (intensity >= GPntile[i])) {
      return(names(GPntile[i]))
    }
  }
}


# Normalize a numeric vector to the sum of its components
normalize <- function (vec) {
  newVec = vec/sum(vec)
  return(newVec)
}


# Generate a table counting up instances of each annotated glycan composition in an annotated PSM table
# Then generate a table counting up how many of each glycan composition are in each quantile group
# Then clean up the output table and generate normalized glycan composition frequencies for each quantile group
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
