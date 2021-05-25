#!/usr/bin/Rscript
library(eQTLpipeline)
library(seqQTL)
library(preprocessCore)
library(data.table)
source("~/work/HS_heart/R/functions/vif_function.R")
source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/eqtl_trans_kinship.R")
source("~/work/new_data_RI_translational_efficiency/Franzi/R/permutations_SDP.R")


res_dir="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/permutations/"

setwd("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Result_tables/20180213_Lv_polyA/")

load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/Lv_polyA_tables.RData")
load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Input/Lv_polyA/Lv_polyA_input.RData")

load(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/filtered_genesets_ribotaper_fpkm.RData")

expr = pheno[which(rownames(pheno) %in% genesset_lv),]
positions = positions[which(positions$gene_id %in% genesset_lv),]
actual.eqtls2 = eQTL_cis

esnps_cis = find.eSNPs_new(expr, covariates, positions, snp.pos, genotypes.file, dir=paste(res_dir, "Lv_polyA_cis", sep=""), threshold=1, min.perm=1000, max.perm=10000, exit.criterion=15, actual.eqtls=actual.eqtls2, only_cis=F)

write.table(esnps_cis, file = paste(res_dir, "Lv_polyA_esnps_cis.txt", sep=""), sep="\t", quote=F, row.names=F)

##########


actual.eqtls3 = eQTL_trans
esnps_trans = find.eSNPs_new(expr, covariates, positions, snp.pos, genotypes.file, dir=paste(res_dir, "Lv_polyA_trans", sep=""), threshold=1, min.perm=1000, max.perm=10000, exit.criterion=15, actual.eqtls=actual.eqtls3, only_cis=F)

write.table(esnps_trans, file = paste(res_dir, "Lv_polyA_esnps_trans.txt", sep=""), sep="\t", quote=F, row.names=F)



save(esnps_cis, file = paste(res_dir, "Lv_polyA_esnps_cis.RData", sep=""))save(esnps_trans, file = paste(res_dir, "Lv_polyA_esnps_trans.RData", sep=""))

