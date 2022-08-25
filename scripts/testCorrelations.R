#!/usr/bin/env Rscript

rm(list=ls())
infileTree <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_helpers/megaphyloTreeSelectedSpecies.outgroup.tre"
infileCounts <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/allCountsPerSpecies.csv"
infileCountsForTree <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/allCountsPerSpeciesSubsetByTree.csv"
outfilePrefix <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/correlationAnalysis"
rm(list=ls())


## extract the path to the script itself
myarg <- commandArgs(trailingOnly = FALSE)
scriptName <- gsub("--file=", "", grep("--file=", myarg, value = TRUE))
scriptDir <- dirname(scriptName)

## arguments from commandline
argPos <- grep("--args", myarg, fixed = TRUE)
infileTree <- as.character(myarg[argPos+1])
infileCounts <- as.character(myarg[argPos+2])
infileCountsForTree <- as.character(myarg[argPos+3])
outfilePrefix <- as.character(myarg[argPos+4])

## load packages silently
suppressPackageStartupMessages({
  library("phytools")
  library("ape")
  library("vegan")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
})



myTree <- read.tree(infileTree)
myDataForTree <- read.csv(infileCountsForTree, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
myData <- read.csv(infileCounts, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
colSums(myData)

####################################################################
# Normalize data by number of searched genes
dataCols <- setdiff(colnames(myData), c("SEARCHED", "SEARCHED150"))
normWith150 <- grep("150$", dataCols, value = TRUE)
normWith250 <- grep("150$", dataCols, value = TRUE, invert = TRUE)

for (curCol in normWith150) {
  myDataForTree[[curCol]] <- myDataForTree[[curCol]]/myDataForTree$SEARCHED150
  myData[[curCol]] <- myData[[curCol]]/myData$SEARCHED150
}

for (curCol in normWith250) {
  myDataForTree[[curCol]] <- myDataForTree[[curCol]]/myDataForTree$SEARCHED
  myData[[curCol]] <- myData[[curCol]]/myData$SEARCHED
}

myDataForTree <- myDataForTree[,dataCols]
myData <- myData[,dataCols]

####################################################################
# Get distance matrices
myDataDistances <- list()
myDataForTreeDistances <- list()
phyloDist <- cophenetic.phylo(myTree)
phyloDist <- phyloDist[rownames(myDataForTree), rownames(myDataForTree)]

for (curCol in dataCols) {
  myDataDistances[[curCol]] <- as.matrix(dist(myData[[curCol]]))
  colnames(myDataDistances[[curCol]]) <- rownames(myData)
  rownames(myDataDistances[[curCol]]) <- rownames(myData)
  myDataForTreeDistances[[curCol]] <- as.matrix(dist(myDataForTree[[curCol]]))
  colnames(myDataForTreeDistances[[curCol]]) <- rownames(myDataForTree)
  rownames(myDataForTreeDistances[[curCol]]) <- rownames(myDataForTree)
}

pValuesMantel <- matrix(NA, nrow = ncol(myData), ncol = ncol(myData), dimnames = list(colnames(myData), colnames(myData)))
pValuesMantelForTree <- matrix(NA, nrow = ncol(myDataForTree)+1, ncol = ncol(myDataForTree)+1, dimnames = list(c(colnames(myDataForTree), "phylogeny"), c(colnames(myDataForTree), "phylogeny")))
pValuesMantelForTreeCorrected <- matrix(NA, nrow = ncol(myDataForTree), ncol = ncol(myDataForTree), dimnames = list(colnames(myDataForTree), colnames(myDataForTree)))

rValuesMantel <- matrix(NA, nrow = ncol(myData), ncol = ncol(myData), dimnames = list(colnames(myData), colnames(myData)))
rValuesMantelForTree <- matrix(NA, nrow = ncol(myDataForTree)+1, ncol = ncol(myDataForTree)+1, dimnames = list(c(colnames(myDataForTree), "phylogeny"), c(colnames(myDataForTree), "phylogeny")))
rValuesMantelForTreeCorrected <- matrix(NA, nrow = ncol(myDataForTree), ncol = ncol(myDataForTree), dimnames = list(colnames(myDataForTree), colnames(myDataForTree)))

####################################################################
# Mantel test

nPerm <- 10000
for (i in 1:(length(dataCols)-1)) {
  first <- dataCols[i]
  f.print.message("Comparing", first, "with all others")
  for (j in (i+1):length(dataCols)) {
    second <- dataCols[j]
    #f.print.message("Comparing", first, "with", second, "all data")
    temp <- mantel(myDataDistances[[first]], myDataDistances[[second]], permutations = nPerm)
    pValuesMantel[first, second] <- temp$signif
    rValuesMantel[first, second] <- temp$statistic
    #f.print.message("Comparing", first, "with", second, "subset data")
    temp <- mantel(myDataForTreeDistances[[first]], myDataForTreeDistances[[second]], permutations = nPerm)
    pValuesMantelForTree[first, second] <- temp$signif
    rValuesMantelForTree[first, second] <- temp$statistic
    #f.print.message("Comparing", first, "with", second, "subset data, partial mantel")
    temp <- mantel.partial(myDataForTreeDistances[[first]], myDataForTreeDistances[[second]], phyloDist, method="pearson", permutations = nPerm)
    pValuesMantelForTreeCorrected[first, second] <- temp$signif
    rValuesMantelForTreeCorrected[first, second] <- temp$statistic
  }
  temp <- mantel(myDataForTreeDistances[[first]], phyloDist, permutations = nPerm)
  pValuesMantelForTree[first, "phylogeny"] <- temp$signif
  rValuesMantelForTree[first, "phylogeny"] <- temp$statistic
}
first <- dataCols[length(dataCols)]
temp <- mantel(myDataForTreeDistances[[first]], phyloDist, permutations = nPerm)
pValuesMantelForTree[first, "phylogeny"] <- temp$signif
rValuesMantelForTree[first, "phylogeny"] <- temp$statistic

fdrValuesMantel <- matrix(p.adjust(as.vector(pValuesMantel), "fdr"), ncol = ncol(pValuesMantel))
fdrValuesMantelForTree <- matrix(p.adjust(as.vector(pValuesMantelForTree), "fdr"), ncol = ncol(pValuesMantelForTree))
fdrValuesMantelForTreeCorrected <- matrix(p.adjust(as.vector(pValuesMantelForTreeCorrected), "fdr"), ncol = ncol(pValuesMantelForTreeCorrected))
dimnames(fdrValuesMantel) <- dimnames(pValuesMantel)
dimnames(fdrValuesMantelForTree) <- dimnames(pValuesMantelForTree)
dimnames(fdrValuesMantelForTreeCorrected) <- dimnames(pValuesMantelForTreeCorrected)

write.csv(pValuesMantel, paste0(outfilePrefix, ".mantel.allData.csv"))
write.csv(pValuesMantelForTree, paste0(outfilePrefix, ".mantel.treeData.csv"))
write.csv(pValuesMantelForTreeCorrected, paste0(outfilePrefix, ".mantel.treeData.phyloAccounted.csv"))
write.csv(rValuesMantel, paste0(outfilePrefix, ".mantel.r.allData.csv"))
write.csv(rValuesMantelForTree, paste0(outfilePrefix, ".mantel.r.treeData.csv"))
write.csv(rValuesMantelForTreeCorrected, paste0(outfilePrefix, ".mantel.treeData.r.phyloAccounted.csv"))
write.csv(fdrValuesMantel, paste0(outfilePrefix, ".mantel.fdr.allData.csv"))
write.csv(fdrValuesMantelForTree, paste0(outfilePrefix, ".mantel.fdr.treeData.csv"))
write.csv(fdrValuesMantelForTreeCorrected, paste0(outfilePrefix, ".mantel.treeData.fdr.phyloAccounted.csv"))

#pValuesMantel <- as.matrix(read.csv(paste0(outfilePrefix, ".mantel.allData.csv"), header = TRUE, row.names = 1))
#pValuesMantelForTree <- as.matrix(read.csv(paste0(outfilePrefix, ".mantel.treeData.csv"), header = TRUE, row.names = 1))
#pValuesMantelForTreeCorrected <- as.matrix(read.csv(paste0(outfilePrefix, ".mantel.treeData.phyloAccounted.csv"), header = TRUE, row.names = 1))

pdf(paste0(outfilePrefix, ".mantel.pValue.plot.pdf"), height = 5, width = 5)
plot(na.omit(as.vector(pValuesMantelForTree[rownames(pValuesMantelForTreeCorrected), colnames(pValuesMantelForTreeCorrected)])), na.omit(as.vector(pValuesMantelForTreeCorrected)), xlab = "P-value uncorrected", ylab = "P-value corrected for phylogeny", main = "P-values (partial) mantel tests")
invisible(dev.off())
#cor(na.omit(as.vector(pValuesMantelForTree[rownames(pValuesMantelForTreeCorrected), colnames(pValuesMantelForTreeCorrected)])), na.omit(as.vector(pValuesMantelForTreeCorrected)))

####################################################################
# the fancy test based on: http://blog.phytools.org/2017/08/pearson-correlation-with-phylogenetic.html
### other ---------------------------------
## correlation between x & y
# r.xy <- cov2cor(obj$R)["LRR_I","LRR_II"]
# ## t-statistic & P-value
# t.xy <- r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
# P.xy <- 2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
# P.xy

pValuesLm <- matrix(NA, nrow = ncol(myData), ncol = ncol(myData), dimnames = list(colnames(myData), colnames(myData)))
pValuesLmForTree <- matrix(NA, nrow = ncol(myDataForTree), ncol = ncol(myDataForTree), dimnames = list(colnames(myDataForTree), colnames(myDataForTree)))
pValuesLmForTreeCorrected <- matrix(NA, nrow = ncol(myDataForTree), ncol = ncol(myDataForTree), dimnames = list(colnames(myDataForTree), colnames(myDataForTree)))

for (i in 1:(length(dataCols)-1)) {
  first <- dataCols[i]
  for (j in (i+1):length(dataCols)) {
    second <- dataCols[j]
    mod <- anova(lm(as.formula(paste0(first, "~", second)), data = myData))
    pValuesLm[first, second] <- mod[1, "Pr(>F)"]
    mod <- anova(lm(as.formula(paste0(first, "~", second)), data = myDataForTree))
    pValuesLmForTree[first, second] <- mod[1, "Pr(>F)"]
    species <- rownames(myDataForTree)
    fit <- nlme::gls(as.formula(paste0(first, "~", second)), data = myDataForTree, correlation = corBrownian(1, myTree, ~ species))
    mod <- anova(fit)
    pValuesLmForTreeCorrected[first, second] <- mod[nrow(mod), "p-value"]
  }
}

write.csv(pValuesLm, paste0(outfilePrefix, ".linearModel.allData.csv"))
write.csv(pValuesLmForTree, paste0(outfilePrefix, ".linearModel.treeData.csv"))
write.csv(pValuesLmForTreeCorrected, paste0(outfilePrefix, ".linearModel.treeData.phyloAccounted.csv"))

pdf(paste0(outfilePrefix, ".linearModel.pValue.plot.pdf"), height = 5, width = 5)
plot(na.omit(as.vector(pValuesLmForTree)), na.omit(as.vector(pValuesLmForTreeCorrected)), xlab = "P-value uncorrected", ylab = "P-value corrected for phylogeny", main = "P-values (partial) mantel tests")
invisible(dev.off())
#cor(na.omit(as.vector(pValuesLmForTree)), na.omit(as.vector(pValuesLmForTreeCorrected)))





