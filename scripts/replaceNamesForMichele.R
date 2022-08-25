#!/usr/bin/env Rscript

rm(list=ls())
infileSpecies <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_helpers/speciesInMegaphyloTree.txt"
infileCountsAtBeginning <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/numberOfGenesAtSearchBegin.csv"
infileCountsAtBeginning150 <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/numberOfGenesAtSearchBegin_150.csv"
infileLRRcounts <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/proteinsPerSpecies.csv"
infileNLRcounts <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/lrrAndNbArc_candidates_per_species.csv"
infileRLPcounts <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/LRR_RLP_candidates_per_species.csv"
infileRLP150counts <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/LRR_RLP_candidates_per_species_150.csv"
infileNbArccounts <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/onlyNbArc_candidates_per_species.csv"
infileNbArc150counts <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/onlyNbArc_candidates_150_per_species.csv"
infileLysRLKcount <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/LysM_RLK_candidates_per_species.csv"
infileLysRLPcount <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/LysM_RLP_candidates_per_species.csv"
infileLysRLK150count <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/LysM_RLK_candidates_150_per_species.csv"
infileLysRLP150count <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/LysM_RLP_candidates_150_per_species.csv"
outfileNameNoSubset <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/allCountsPerSpecies.csv"
outfileName <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/allCountsPerSpeciesSubsetByTree.csv"
rm(list=ls())


## extract the path to the script itself
myarg <- commandArgs(trailingOnly = FALSE)
scriptName <- gsub("--file=", "", grep("--file=", myarg, value = TRUE))
scriptDir <- dirname(scriptName)

## arguments from commandline
argPos <- grep("--args", myarg, fixed = TRUE)
infileSpecies <- as.character(myarg[argPos+1])
infileCountsAtBeginning <- as.character(myarg[argPos+2])
infileCountsAtBeginning150 <- as.character(myarg[argPos+3])
infileLRRcounts <- as.character(myarg[argPos+4])
infileNLRcounts <- as.character(myarg[argPos+5])
infileRLPcounts <- as.character(myarg[argPos+6])
infileRLP150counts <- as.character(myarg[argPos+7])
infileNbArccounts <- as.character(myarg[argPos+8])
infileNbArc150counts <- as.character(myarg[argPos+9])
infileLysRLKcount <- as.character(myarg[argPos+10])
infileLysRLPcount <- as.character(myarg[argPos+11])
infileLysRLK150count <- as.character(myarg[argPos+12])
infileLysRLP150count <- as.character(myarg[argPos+13])
outfileNameNoSubset <- as.character(myarg[argPos+14])
outfileName <- as.character(myarg[argPos+15])

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
indCounts <- list(
  beg = read.csv(infileCountsAtBeginning, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
  beg150 = read.csv(infileCountsAtBeginning150, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
  lrr = read.csv(infileLRRcounts, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
  nlr = read.csv(infileNLRcounts, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
  rlp = read.csv(infileRLPcounts, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
  rlp150 = read.csv(infileRLP150counts, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
  nbArc = read.csv(infileNbArccounts, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
  nbArc150 = read.csv(infileNbArc150counts, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
  lysRlk = read.csv(infileLysRLKcount, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
  lysRlp = read.csv(infileLysRLPcount, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
  lysRlk150 = read.csv(infileLysRLK150count, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
  lysRlp150 = read.csv(infileLysRLP150count, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
)

colnames(indCounts$beg) <- "SEARCHED"
colnames(indCounts$beg150) <- "SEARCHED150"
colnames(indCounts$nlr) <- "NLR"
colnames(indCounts$nbArc) <- "NbArc"
colnames(indCounts$nbArc150) <- "NbArc150"
colnames(indCounts$lysRlk) <- "LysMrlk"
colnames(indCounts$lysRlp) <- "LysMrlp"
colnames(indCounts$lysRlk150) <- "LysMrlk150"
colnames(indCounts$lysRlp150) <- "LysMrlp150"
colnames(indCounts$rlp)[colnames(indCounts$rlp)=="count"] <- "anyRLP"
colnames(indCounts$rlp150)[colnames(indCounts$rlp150)=="count"] <- "anyRLP"
colnames(indCounts$rlp150) <- paste0(colnames(indCounts$rlp150), "150")
allSpecies <- Reduce(union, lapply(indCounts, rownames))
allColumns <- Reduce(union, lapply(indCounts, colnames))
allCounts <- matrix(0, nrow = length(allSpecies), ncol = length(allColumns), dimnames = list(allSpecies, allColumns))
for (toAdd in names(indCounts)) {
  temp <- indCounts[[toAdd]]
  for (curCol in colnames(temp)) {
    allCounts[rownames(temp), curCol] <- temp[,curCol]
  }
}

write.csv(allCounts, outfileNameNoSubset)

################################################################################################
### match names
originalSpeciesShort <- sapply(originalSpecies, function(x) {temp <- unlist(strsplit(x, '_', fixed = TRUE)); return(paste0(substr(temp[1],1,1), temp[2]))})
originalSpeciesShort[originalSpeciesShort == "Pintegrifolia"] <- "Pinflata"
originalSpeciesShort[originalSpeciesShort == "Amajus"] <- "AmajusL"
originalSpeciesShort[originalSpeciesShort == "Ctrifoliata"] <- "Ptrifoliata"
originalSpeciesShort[originalSpeciesShort == "Abuxifolia"] <- "Abuxfoliata"
missing <- setdiff(originalSpeciesShort, rownames(allCounts))
if (length(missing) > 0) {
  f.print.message("ERROR: not everything here, exiting.")
  quit("no", 0)
}
names(originalSpecies) <- originalSpeciesShort
allCounts <- allCounts[names(originalSpecies),]
rownames(allCounts) <- originalSpecies
write.csv(allCounts, outfileName)

