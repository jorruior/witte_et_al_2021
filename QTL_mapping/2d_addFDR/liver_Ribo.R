#!/usr/bin/Rscript

library(GenomicRanges)
library(eQTLpipeline)
library(data.table)
source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/SeqQTL_oldCluster.R")
dir="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Result_tables/"


load(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Result_tables/local_sdps.RData")
load(paste(dir, "20180213_Liver_Ribo/image.RData", sep=""))
load(file=paste(dir, "20180213_Liver_Ribo/eqtls_sorted.RData", sep=""))
eQTL_results_Liver_Ribo  = data.table(eQTL_results_liver_ribo)
setkey(eQTL_results_Liver_Ribo ,gene,chrom,snp_pos)

unknowns = NULL
for(counter in c(1:22))
{

	load(file =  paste(dir, "20180213_Liver_Ribo/unknowns_", counter, ".RData", sep=""))
	unknowns = rbind(unknowns, t(locals_unknowns))
}



test[ids_unknown,] = unknowns
kind=rep("NA", nrow(eQTL_results_Liver_Ribo))
setkey(eQTL_results_Liver_Ribo,gene,chrom,snp_pos)

kind[eQTL_results_Liver_Ribo[[9]] != test[,1]] = "trans"
kind[eQTL_results_Liver_Ribo[[9]] == test[,1]] = "cis"
kind[eQTL_results_Liver_Ribo[[1]] == eQTL_results_Liver_Ribo[[14]]] = "cis"

eQTL_results_Liver_Ribo$kind = kind
eQTL_results_Liver_Ribo$chr_sdp_test = test[,1]
eQTL_results_Liver_Ribo$sdp_test = test[,2]
setkey(eQTL_results_Liver_Ribo, kind)
eQTL_trans = eQTL_results_Liver_Ribo[kind == "trans"]
eQTL_cis = eQTL_results_Liver_Ribo[kind == "cis"]

load(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/filtered_genesets_ribotaper_fpkm.RData")


eQTL_cis = eQTL_cis[gene %in% genesset_li]
eQTL_trans = eQTL_trans[gene %in% genesset_li]



eQTL_trans$FDR_new = p.adjust(eQTL_trans[[5]], method="BH")
eQTL_cis$FDR_new = p.adjust(eQTL_cis[[5]], method="BH")

setkey(eQTL_trans, FDR_new)
setkey(eQTL_cis, FDR_new)

eQTL_trans.sign = eQTL_trans[FDR_new <= 0.1]
#12 pairs
#10 unique genes
setkey(eQTL_trans.sign, gene)

write.table(eQTL_trans.sign, file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/distant_eQTLs_Liver_Ribo_sign.txt", sep="\t", col.names=T, row.names=F, quote=F)

eQTL_cis.sign = eQTL_cis[FDR_new <= 0.05]
#6541 pairs
#1268 unique genes
setkey(eQTL_cis.sign, gene)
write.table(eQTL_cis.sign, file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/local_eQTLs_Liver_Ribo_sign.txt", sep="\t", col.names=T, row.names=F, quote=F)

save.image("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/Liver_Ribo_tables.RData")

