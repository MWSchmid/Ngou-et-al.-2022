#!/usr/bin/env Rscript

rm(list=ls())
infileSpecies <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_helpers/speciesForTree.txt"
infileTree <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_helpers/megaphylo/tree/PhytoPhylo.tre"
outfileName <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_helpers/speciesInMegaphyloTree.txt"
outfileNameMissing <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_helpers/speciesMissingInMegaphyloTree.txt"
rm(list=ls())

## extract the path to the script itself
myarg <- commandArgs(trailingOnly = FALSE)
scriptName <- gsub("--file=", "", grep("--file=", myarg, value = TRUE))
scriptDir <- dirname(scriptName)

## arguments from commandline
argPos <- grep("--args", myarg, fixed = TRUE)
infileSpecies <- as.character(myarg[argPos+1])
infileTree <- as.character(myarg[argPos+2])
outfileName <- as.character(myarg[argPos+3])
outfileNameMissing <- as.character(myarg[argPos+4])

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
originalSpecies <- scan(infileSpecies, what = "character", sep = '\n')
originalSpecies <- gsub(" C.*", "", originalSpecies)
originalSpeciesWithUnderscore <- gsub(" ", "_", originalSpecies)
myTree <- read.tree(infileTree)
treeOrder <- myTree$tip.label
presentSpecies <- intersect(originalSpeciesWithUnderscore, treeOrder)
f.print.message("found", length(presentSpecies), "species.")
write(presentSpecies, file = outfileName)
missingSpecies <- setdiff(treeOrder, presentSpecies)
write(missingSpecies, file = outfileNameMissing)