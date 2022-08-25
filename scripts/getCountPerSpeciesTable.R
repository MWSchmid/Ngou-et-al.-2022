#!/usr/bin/env Rscript

rm(list=ls())
ncbiInput <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_helpers/allSpeciesList.txt"
treeFile <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_helpers/Revised_tree.phy"
treeFile <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_helpers/Revised_tree_marc.phy"
infileName <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/candsToArabidopsisSubgroup.csv"
outfileName <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/proteinsPerSpecies.csv"
rdKinasesInfile <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/RD_kinases.txt"
rm(list=ls())

## extract the path to the script itself
myarg <- commandArgs(trailingOnly = FALSE)
scriptName <- gsub("--file=", "", grep("--file=", myarg, value = TRUE))
scriptDir <- dirname(scriptName)

## arguments from commandline
argPos <- grep("--args", myarg, fixed = TRUE)
ncbiInput <- as.character(myarg[argPos+1])
treeFile <- as.character(myarg[argPos+2])
infileName <- as.character(myarg[argPos+3])
outfileName <- as.character(myarg[argPos+4])

getRDkinases <- FALSE
if ("--RDkinases" %in% myarg){
  getRDkinases <- TRUE
  rdKinasesInfile <- myarg[which(myarg=="--RDkinases")+1]
}

## load packages silently
suppressPackageStartupMessages({
  library("methods")
  library("ape")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
})

################################################################################################
### load data
cat("reading...\n")
myData <- read.csv(infileName, header = TRUE, stringsAsFactors = FALSE)
myData$species <- sapply(myData$candidateNameOnly, function(x) unlist(strsplit(x, '_'))[1])
originalSpecies <- scan(ncbiInput, what = "character", sep = '\n')
wrongNames <- setdiff(myData$species, originalSpecies)
if (length(wrongNames) > 0) {
  cat("ERROR with species names!\n")
  quit("no", 0)
}

myTree <- read.tree(treeFile)
treeOrder <- myTree$tip.label
length(treeOrder)
length(originalSpecies)
#cat(paste0(gsub(".", " ", treeOrder, fixed = TRUE), collapse = '\n'), '\n')
treeOrder <- gsub("Cannabis.sativa", "C.annabissativa", treeOrder, fixed = TRUE)
treeOrder <- gsub("Ceratodon.purpureus", "C.eratodonpurpureus", treeOrder, fixed = TRUE)
treeOrder <- gsub("Citrus.sinensis", "C.itrussinensis", treeOrder, fixed = TRUE)
treeOrder <- gsub("Pyrus.x.bretschneideri", "P.bretschneideri", treeOrder, fixed = TRUE)
treeOrder <- gsub("Atalantia.buxifolia", "Atalantia.buxfoliata", treeOrder, fixed = TRUE)
treeOrder <- gsub("Saccharum.hybrid.cultivar", "Saccharum.spontaneumhybrid", treeOrder, fixed = TRUE)
treeOrder <- gsub("Dioscorea.cayenensis.subsp..rotundata", "Dioscorea.rotundata", treeOrder, fixed = TRUE)
treeOrder <- gsub("Dioscorea.cayenensis", "Dioscorea.rotundata", treeOrder, fixed = TRUE)
treeOrder <- gsub("Antirrhinum.majus", "Antirrhinum.majusL", treeOrder, fixed = TRUE)
treeOrder <- gsub("Dorcoceras.hygrometricum", "Boea.hygrometrica", treeOrder, fixed = TRUE)
treeOrder <- gsub("Erythranthe.guttata", "Mimulus.guttatus", treeOrder, fixed = TRUE)
treeOrder <- gsub("Meniocus.linifolius", "Alyssum.linifolium", treeOrder, fixed = TRUE)
treeOrder <- gsub("Citrus.trifoliata", "Poncirus.trifoliata", treeOrder, fixed = TRUE)
treeOrder <- gsub("Petunia.integrifolia.subsp..inflata", "Petunia.inflata", treeOrder, fixed = TRUE)
treeOrder <- gsub("Petunia.integrifolia", "Petunia.inflata", treeOrder, fixed = TRUE)
treeOrder <- gsub("Citrus.cavaleriei", "Citrus.ichangensis", treeOrder, fixed = TRUE)
treeOrder <- gsub("Cinnamomum.micranthum.f..kanehirae", "Cinnamomum.kanehirae", treeOrder, fixed = TRUE) # the C. micranthum is without the f. kanehirae, we deleted the plain C. micranthum because one couldn't not include both in the tree
treeOrder <- gsub("Citrus.maxima", "Citrus.grandis", treeOrder, fixed = TRUE)
treeOrder <- gsub("Citrus.hindsii", "Fortunella.hindsii", treeOrder, fixed = TRUE)
treeOrder <- gsub(".x.", ".x", treeOrder, fixed = TRUE)
treeOrder <- sapply(treeOrder, function(x) {temp <- unlist(strsplit(x, '.', fixed = TRUE)); return(paste0(substr(temp[1],1,1), temp[2]))})
#head(sort(table(treeOrder), decreasing = TRUE))
myTree$tip.label <- sapply(treeOrder, function(x) paste0(substr(x, 1, 1), ". ", substr(x, 2, nchar(x))))
missing <- union(setdiff(treeOrder, originalSpecies), setdiff(originalSpecies, treeOrder))
if (length(missing)>0) {
  cat("ERROR!\n")
  quit("no", 0)
}

################################################################################################
# check what's additional in the tree
count <- rep(1, nrow(myData))
countPerSpeciesAndGroupLong <- aggregate(count, by = list(species = myData$species, subgroup = myData$subgroup), sum)
countPerSpeciesAndGroup <- f.long.to.wide(countPerSpeciesAndGroupLong, "species", "subgroup", "x")
missingSpecies <- setdiff(treeOrder, rownames(countPerSpeciesAndGroup))
toAdd <- matrix(0, nrow = length(missingSpecies), ncol = ncol(countPerSpeciesAndGroup), dimnames = list(missingSpecies, colnames(countPerSpeciesAndGroup)))
countPerSpeciesAndGroup <- rbind(countPerSpeciesAndGroup, toAdd)
countPerSpeciesAndGroup[is.na(countPerSpeciesAndGroup)] <- 0
colnames(countPerSpeciesAndGroup) <- gsub("-", "_", colnames(countPerSpeciesAndGroup))
forOrder <- c(
  "LRR_I",
  "LRR_II",
  "LRR_III",
  "LRR_IV",
  "LRR_V",
  "LRR_VI_1",
  "LRR_VI_2",
  "LRR_VII",
  "LRR_VIII_1",
  "LRR_VIII_2",
  "LRR_IX",
  "LRR_Xa",
  "LRR_Xb",
  "LRR_XI",
  "LRR_XII",
  "LRR_XIIIa",
  "LRR_XIIIb",
  "LRR_XIV",
  "LRR_XV",
  "LRR_XVI"
)
countPerSpeciesAndGroup <- countPerSpeciesAndGroup[treeOrder, forOrder]
write.csv(countPerSpeciesAndGroup, outfileName)
write.tree(myTree, gsub(".csv", ".tree", outfileName, fixed = TRUE))

################################################################################################
# classify into RD kinases
if (!getRDkinases) {
  cat("no RD kinases provided. Exiting...\n")
  quit("no", 0)
}

rdKinases <- scan(rdKinasesInfile, what = "character")
myData$kinaseType <- "nonRD"
myData$kinaseType[myData$candidate %in% rdKinases] <- "RD"
table(myData$kinaseType)
for (curType in unique(myData$kinaseType)) {
  temp <- subset(myData, kinaseType == curType)
  count <- rep(1, nrow(temp))
  countPerSpeciesAndGroupLong <- aggregate(count, by = list(species = temp$species, subgroup = temp$subgroup), sum)
  countPerSpeciesAndGroup <- f.long.to.wide(countPerSpeciesAndGroupLong, "species", "subgroup", "x")
  missingSpecies <- setdiff(treeOrder, rownames(countPerSpeciesAndGroup))
  toAdd <- matrix(0, nrow = length(missingSpecies), ncol = ncol(countPerSpeciesAndGroup), dimnames = list(missingSpecies, colnames(countPerSpeciesAndGroup)))
  countPerSpeciesAndGroup <- rbind(countPerSpeciesAndGroup, toAdd)
  countPerSpeciesAndGroup[is.na(countPerSpeciesAndGroup)] <- 0
  colnames(countPerSpeciesAndGroup) <- gsub("-", "_", colnames(countPerSpeciesAndGroup))
  countPerSpeciesAndGroup <- countPerSpeciesAndGroup[treeOrder, forOrder]
  write.csv(countPerSpeciesAndGroup, gsub("\\.csv$", paste0(".", curType, ".csv"), outfileName))
}




