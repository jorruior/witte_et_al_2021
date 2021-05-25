#Script to estimate heritability based on the variance of replicates (BXH13 and BXH12). Power was calculated using this online tool: https://github.com/Dashbrook/BXD_power_calculator_app
library(DESeq2)

rep = read.table("replicate_counts.txt") #Table of counts from the selected replicates (3 vs 3)

va <- function(x,y){
  v = var(c(x,y))
  return(v)
}

ve <- function(x,y){
  v = (var(x) + var(y)) / 2
  return(v)
}

her <- function(x){
  h = (0.5*x["va"])/((0.5*x["va"])+x["ve"])
  return(h)
}

dds_polyA=DESeqDataSetFromMatrix(
  countData = rep,
  colData = data.frame(sample_id = colnames(rep), condition = c(1,1,1,2,2,2)),
  design = ~ condition)

dds <- DESeq(dds_polyA)
dds <- estimateSizeFactors(dds)
ndds <- data.frame(counts(dds, normalized=TRUE))

ndds["mean"] = apply(ndds,1,mean)
ndds = ndds[ndds$mean > 50,]
ndds["va"] = apply(ndds,1,function(x) va(x[c("BXH12_2","BXH12_3","BXH12_4")],x[c("BXH13_2","BXH13_3","BXH13_4")]))
ndds["ve"] = apply(ndds,1,function(x) ve(x[c("BXH12_2","BXH12_3","BXH12_4")],x[c("BXH13_2","BXH13_3","BXH13_4")]))
ndds["v12"] = apply(ndds[,c("BXH12_2","BXH12_3","BXH12_4")],1,cv)
ndds["v13"] = apply(ndds[,c("BXH13_2","BXH13_3","BXH13_4")],1,cv)
ndds["v1213"] = apply(ndds[,c("BXH13_2","BXH13_3","BXH13_4","BXH12_2","BXH12_3","BXH12_4")],1,cv)
ndds["h2t"] = apply(ndds,1,function(x) her(x))
