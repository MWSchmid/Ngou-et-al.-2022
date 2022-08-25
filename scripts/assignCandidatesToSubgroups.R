#!/usr/bin/env Rscript

rm(list=ls())
infileName <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/candsToArabidopsisBlast.tsv.gz"
outfileName <- "/media/mwschmid/myData/MWSchmid/JJ_Bruno/GitIgnore_results/candsToArabidopsisSubgroup.csv"
rm(list=ls())

## extract the path to the script itself
myarg <- commandArgs(trailingOnly = FALSE)
scriptName <- gsub("--file=", "", grep("--file=", myarg, value = TRUE))
scriptDir <- dirname(scriptName)

## arguments from commandline
argPos <- grep("--args", myarg, fixed = TRUE)
infileName <- as.character(myarg[argPos+1])
outfileName <- as.character(myarg[argPos+2])

## load packages silently
suppressPackageStartupMessages({
  library("parallel")
})

################################################################################################
### load data
cat("reading...\n")
myData <- read.table(gzfile(infileName), header = FALSE, stringsAsFactors = FALSE, sep = '\t', col.names = c("query", "subject", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "score"))
myData <- myData[,c("query", "subject", "score")]
#myData$subgroup <- sapply(myData$subject, function(x) unlist(strsplit(x, '__'))[1])

# split by query and then use mclapply(), index split would be faster
cat("splitting...\n")
splitIndices <- split(1:nrow(myData), myData$query)

# select the subgroup based on the highest average bitscore
f.select.subgroup <- function(x) {
  query <- x$query[1]
  # only best hit per pair, if there are more than one
  subjectCount <- table(x$subject)
  if (sum(subjectCount > 1) > 0) {
    temp <- split(x, x$subject)
    temp <- lapply(temp, function(y) {if (nrow(y) > 1) {return(y[which.max(y$score),])} else {return (y)}})
    temp <- as.data.frame(do.call("rbind", temp), stringsAsFactors = FALSE)
    temp$score <- as.numeric(temp$score)
  } else {
    temp <- x
  }
  temp$subgroup <- sapply(temp$subject, function(x) unlist(strsplit(x, '__'))[1]) # move it here do profit from parallel
  meanPerGroup <- aggregate(temp$score, by = list(group = temp$subgroup), mean)
  selectedGroup <- meanPerGroup$group[which.max(meanPerGroup$x)]
  out <- c(query, selectedGroup)
  return(out)
}

cat("extracting...\n")
results <- mclapply(splitIndices, function(x) f.select.subgroup(myData[x,]), mc.cores = 30)
results <- as.data.frame(do.call("rbind", results), stringsAsFactors = FALSE)
colnames(results) <- c("candidate", "subgroup")
results$subgroup <- gsub("at_sg_", "", results$subgroup)
results$candidateNameOnly <- sapply(results$candidate, function(x) unlist(strsplit(x, ":"))[1])
table(results$subgroup)
write.csv(results, outfileName, quote = FALSE, row.names = FALSE)


