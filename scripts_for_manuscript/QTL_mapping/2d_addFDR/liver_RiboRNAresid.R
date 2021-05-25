#!/usr/bin/Rscript

library(GenomicRanges)
library(eQTLpipeline)
library(data.table)
source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/SeqQTL_oldCluster.R")
dir="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Result_tables/"


load(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Result_tables/local_sdps.RData")
load(paste(dir, "20180403_Liver_Ribo_RNA_resid/image.RData", sep=""))
load(file=paste(dir, "20180403_Liver_Ribo_RNA_resid/eqtls_sorted.RData", sep=""))
eQTL_results_Liver_RiboRNAresid  = data.table(eQTL_results_liver_riboRnaResid)
setkey(eQTL_results_Liver_RiboRNAresid ,gene,chrom,snp_pos)

unknowns = NULL
for(counter in c(1:22))
{

	load(file =  paste(dir, "20180403_Liver_Ribo_RNA_resid/unknowns_", counter, ".RData", sep=""))
	unknowns = rbind(unknowns, t(locals_unknowns))
}
test2 = as.matrix(test)
test2[as.numeric(ids_unknown),] = unknowns


kind=rep("NA", nrow(eQTL_results_Liver_RiboRNAresid))
setkey(eQTL_results_Liver_RiboRNAresid,gene,chrom,snp_pos)

kind[eQTL_results_Liver_RiboRNAresid[[9]] != test[,1]] = "trans"
kind[eQTL_results_Liver_RiboRNAresid[[9]] == test[,1]] = "cis"
kind[eQTL_results_Liver_RiboRNAresid[[1]] == eQTL_results_Liver_RiboRNAresid[[14]]] = "cis"

eQTL_results_Liver_RiboRNAresid$kind = kind
eQTL_results_Liver_RiboRNAresid$chr_sdp_test = test[,1]
eQTL_results_Liver_RiboRNAresid$sdp_test = test[,2]
setkey(eQTL_results_Liver_RiboRNAresid, kind)
eQTL_trans = eQTL_results_Liver_RiboRNAresid[kind == "trans"]
eQTL_cis = eQTL_results_Liver_RiboRNAresid[kind == "cis"]

load(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/filtered_genesets_ribotaper_fpkm.RData")


eQTL_cis = eQTL_cis[gene %in% genesset_li]
eQTL_trans = eQTL_trans[gene %in% genesset_li]



eQTL_trans$FDR_new = p.adjust(eQTL_trans[[5]], method="BH")
eQTL_cis$FDR_new = p.adjust(eQTL_cis[[5]], method="BH")

setkey(eQTL_trans, FDR_new)
setkey(eQTL_cis, FDR_new)

eQTL_trans.sign = eQTL_trans[FDR_new <= 0.1]
#2 pairs
#2 unique genes
setkey(eQTL_trans.sign, gene)

write.table(eQTL_trans.sign, file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/distant_eQTLs_Liver_RiboRNAresid_sign.txt", sep="\t", col.names=T, row.names=F, quote=F)

eQTL_cis.sign = eQTL_cis[FDR_new <= 0.05]
#749 pairs
#215 unique genes
setkey(eQTL_cis.sign, gene)
write.table(eQTL_cis.sign, file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/local_eQTLs_Liver_RiboRNAresid_sign.txt", sep="\t", col.names=T, row.names=F, quote=F)

save.image("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/Liver_RiboRNAresid_tables.RData")

