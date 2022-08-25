# sort best kinase
# input: genename, startdomani, enddomain, score
# Wyler M., 

args <- commandArgs(TRUE)
TAB <- args[1]

suppressPackageStartupMessages({
  library("parallel")
  library("data.table")
})
tabella <- as.data.frame(fread(TAB, sep = '\t'))

# 
indexSplit <- split(1:nrow(tabella), tabella$V1)
f.select.hit <- function(x) {
  if (!(is.matrix(x)|is.data.frame(x))) {
    out <- x
  } else {
    out <- x[which.max(x[,4]),]
  }
  return(out)
}
temp <- lapply(indexSplit, function(x) f.select.hit(tabella[x,]))
temp <- do.call("rbind", temp)
write.table(temp[,1:3], col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

quit("no", 0)

# original from Michele, picks the longes, not highest scoring (it^s usually similar though)
nuovaTabella <- as.data.frame(matrix(ncol = 3, nrow = length(unique(tabella$V1))))
nuovaTabella$V1 <- unique(tabella$V1)

for (ID in unique(tabella$V1)){
  subsetTab <- tabella[tabella$V1 == ID, ]
  if (nrow(subsetTab)>1){
    # order by score
    subsetTab <- subsetTab[order(subsetTab$V4, decreasing = T),]
    # pick longest
    subsetTab$LEN <- subsetTab$V3- subsetTab$V2
    write.table(subsetTab[order(subsetTab$LEN, decreasing = T)[1], 1:3],
          col.names = F, 
          quote = F, 
          row.names = F, sep = '\t')
  }else{
    write.table(subsetTab[,1:3],
                col.names = F, 
                quote = F, 
                row.names = F, sep = '\t')
  }
}
