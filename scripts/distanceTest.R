#!/usr/bin/env Rscript

rm(list=ls())
referenceSetFile <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_testing/Slycopersicum_NLRnoLRR.mrnaPos.csv"
testSetFile <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_testing/Slycopersicum_RLP.mrnaPos.csv"
samplingFile <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_testing/Slycopersicum.toIncludeWhileSampling.short.mrnaPos.csv"
outfilePrefix <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_testing/Slycopersicum.RLP.to.NLRnoLRR"
referenceSetFile <- "/media/localData/tempJJB/Ntabacum_NLRnoLRR.mrnaPos.csv"
testSetFile <- "/media/localData/tempJJB/Ntabacum_LRR-III.mrnaPos.csv"
samplingFile <- "/media/localData/tempJJB/Ntabacum.toIncludeWhileSampling.long.mrnaPos.csv"
outfilePrefix <- "/media/localData/tempJJB/test.Ntabacum.LRR-III.to.NLRnoLRR"
rm(list=ls())

## extract the path to the script itself
myarg <- commandArgs(trailingOnly = FALSE)
scriptName <- gsub("--file=", "", grep("--file=", myarg, value = TRUE))
scriptDir <- dirname(scriptName)

## arguments from commandline
argPos <- grep("--args", myarg, fixed = TRUE)
referenceSetFile <- as.character(myarg[argPos+1])
testSetFile <- as.character(myarg[argPos+2])
samplingFile <- as.character(myarg[argPos+3])
outfilePrefix <- as.character(myarg[argPos+4])

## load packages silently
suppressPackageStartupMessages({
  library("methods")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
})

set.seed(333)

################################################################################################
### load data
referenceSet <- read.csv(referenceSetFile, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "pos"))
testSet <- read.csv(testSetFile, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "pos"))
samplingSet <- read.csv(samplingFile, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "pos"))
numTestSet <- nrow(testSet)

################################################################################################
### function for distance calculation
f.get.distance.to.next <- function(dataSetA, dataSetB) { # a is test, b if reference
  commonChroms <- intersect(dataSetA$chrom, dataSetB$chrom)
  numOutsideChromosome <- sum(!(dataSetA$chrom %in% commonChroms))
  minDists <- list()
  for (curChrom in commonChroms) {
    setA <- subset(dataSetA, chrom == curChrom)
    setB <- subset(dataSetB, chrom == curChrom)
    minDists[[curChrom]] <- sapply(setA$pos, function(x) min(abs(x-setB$pos)))
  }
  minDists <- unlist(minDists)
  #out <- c(nrow(dataSetA), numOutsideChromosome, mean(minDists), sd(minDists))
  #out <- c(nrow(dataSetA), numOutsideChromosome, median(minDists), sd(minDists))
  out <- c(nrow(dataSetA), numOutsideChromosome, quantile(minDists, 0.9), sd(minDists))
  return(out)
}

################################################################################################
### get the observed and the samples distances
# keep only chromosomes/scaffolds with at least five entries?
chromCount <- table(samplingSet$chrom)
smallChroms <- names(chromCount)[chromCount < 5]
if (length(smallChroms)) {
  f.print.message("removing", length(smallChroms), "small scaffolds")
  samplingSet <- subset(samplingSet, !(chrom %in% smallChroms))
  referenceSet <- subset(referenceSet, !(chrom %in% smallChroms))
  testSet <- subset(testSet, !(chrom %in% smallChroms))
  numTestSet <- nrow(testSet)
}

observed <- f.get.distance.to.next(testSet, referenceSet); names(observed) <- c("n", "numNotOnChrom", "mean", "sd")
sampled <- list()
for (i in 1:1000) {
  sampledSet <- samplingSet[sample.int(nrow(samplingSet), numTestSet),]
  sampled[[i]] <-f.get.distance.to.next(sampledSet, referenceSet)
}
sampled <- do.call("rbind", sampled)
colnames(sampled) <- names(observed)

################################################################################################
### plot it
allDistances <- sampled[,"mean"]
pValueForTable <- mean(allDistances < observed["mean"], na.rm = TRUE) # with the scaffolds there are some cases that don't share a single scaffold
if (pValueForTable == 0) {
  pValue <- "< 0.001"
} else {
  pValue <- paste0("= ", pValueForTable)
}
pdf(paste0(outfilePrefix, ".density.pdf"), width = 5, height = 5)
forPlot <- density(allDistances, bw = "SJ")
plot(forPlot, main = paste0("Distance to next reference gene, P ", pValue), xlab = "Distance", las = 1)
rug(allDistances, col = f.add.alpha("black", 0.3))
points(observed["mean"], max(forPlot$y), col = "red", pch = 6)
abline(v=observed["mean"], col = "red", lty = 2)
invisible(dev.off())

################################################################################################
### print some metrics
cat(paste(basename(outfilePrefix), observed["n"], round(observed["mean"], 2), round(mean(allDistances), 2), round(observed["mean"]/mean(allDistances), 3), pValueForTable, sep = '\t'), '\n')



