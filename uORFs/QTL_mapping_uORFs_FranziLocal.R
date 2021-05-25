## QTL mapping uORFs

library(eQTLpipeline)
library(data.table)
library(GenomicRanges)
library(grid)

#source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/eqtl_trans_kinship.R")
#source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/SeqQTL_oldCluster.R")




load("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment/20180227_image.RData")
res_dir="/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment/"
sample_names = c("BXH02", "BXH03", "BXH05", "HXB01", "HXB02", "HXB03", "HXB04", "HXB05", "HXB07", "HXB10", "HXB13", "HXB15", "HXB17", "HXB18", "HXB20", "HXB21", "HXB22", "HXB23", "HXB24", "HXB25", "HXB27", "HXB29", "HXB31", "BXH06", "BXH08", "BXH09", "BXH10", "BXH11", "BXH12", "BXH13") 

rownames(length_Norm_ribo_uORF_li) = length_Norm_ribo_uORF_li$ORF_id_tr
colnames(length_Norm_ribo_uORF_li) = gsub("_4", "",gsub("_li", "", colnames(length_Norm_ribo_uORF_li)))
length_Norm_ribo_uORF_li = length_Norm_ribo_uORF_li[,-31]

rownames(length_Norm_ribo_uORF_lv) = length_Norm_ribo_uORF_lv$ORF_id_tr
colnames(length_Norm_ribo_uORF_lv) = gsub("_lv", "", colnames(length_Norm_ribo_uORF_lv))
length_Norm_ribo_uORF_lv = length_Norm_ribo_uORF_lv[,-31]


uORF_ribo_li = length_Norm_ribo_uORF_li[,match(sample_names, colnames(length_Norm_ribo_uORF_li))]
phenotype.file1 = paste(res_dir, "Input_MatrixEQTL/Liver_uORF/phenotypes.txt", sep="")
write.table(uORF_ribo_li, phenotype.file1, sep="\t", quote=F, row.names=T, col.names=T)

uORF_ribo_lv = length_Norm_ribo_uORF_lv[,match(sample_names, colnames(length_Norm_ribo_uORF_lv))]
phenotype.file2 = paste(res_dir, "Input_MatrixEQTL/Lv_uORF/phenotypes.txt", sep="")
write.table(uORF_ribo_lv, phenotype.file2, sep="\t", quote=F, row.names=T, col.names=T)



#Make a uORF position file 
genes_filtered_lv = data.frame(gene_id = all_uORFs_lv[match(rownames(uORF_ribo_lv), all_uORFs_lv$ORF_id_tr),]$ORF_id_tr, chrom = str_split_fixed(all_uORFs_lv[match(rownames(uORF_ribo_lv), all_uORFs_lv$ORF_id_tr),]$ORF_id_gen,"_",3)[,1], start = str_split_fixed(all_uORFs_lv[match(rownames(uORF_ribo_lv), all_uORFs_lv$ORF_id_tr),]$ORF_id_gen, "_",3)[,2], end = str_split_fixed(all_uORFs_lv[match(rownames(uORF_ribo_lv), all_uORFs_lv$ORF_id_tr),]$ORF_id_gen,"_",3)[,3], gene_name=all_uORFs_lv[match(rownames(uORF_ribo_lv), all_uORFs_lv$ORF_id_tr),]$gene_symbol, gene_type=all_uORFs_lv[match(rownames(uORF_ribo_lv), all_uORFs_lv$ORF_id_tr),]$uORF_type)

write.table(genes_filtered_lv, file = paste(res_dir, "Input_MatrixEQTL/Lv_uORF/gene_pos_lv.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)


genes_filtered_li = data.frame(gene_id = all_uORFs_li[match(rownames(uORF_ribo_li), all_uORFs_li$ORF_id_tr),]$ORF_id_tr, chrom = str_split_fixed(all_uORFs_li[match(rownames(uORF_ribo_li), all_uORFs_li$ORF_id_tr),]$ORF_id_gen,"_",3)[,1], start = str_split_fixed(all_uORFs_li[match(rownames(uORF_ribo_li), all_uORFs_li$ORF_id_tr),]$ORF_id_gen, "_",3)[,2], end = str_split_fixed(all_uORFs_li[match(rownames(uORF_ribo_li), all_uORFs_li$ORF_id_tr),]$ORF_id_gen,"_",3)[,3], gene_name=all_uORFs_li[match(rownames(uORF_ribo_li), all_uORFs_li$ORF_id_tr),]$gene_symbol, gene_type=all_uORFs_li[match(rownames(uORF_ribo_li), all_uORFs_li$ORF_id_tr),]$uORF_type)

write.table(genes_filtered_li, file = paste(res_dir, "Input_MatrixEQTL/Liver_uORF/gene_pos_li.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)





#Lv 
pheno = uORF_ribo_lv
phenotype.file = paste(res_dir, "Input_MatrixEQTL/Lv_uORF/phenotypes.txt", sep="")

qtl_dir = "/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/"

genotypes.file =  paste(qtl_dir, "Input/sdp_genotypes.txt", sep="")
covariates.file = NULL
snps_location_file_name = paste(qtl_dir, "Input/snp_pos.txt", sep="")
gene_location_file_name = paste(res_dir, "Input_MatrixEQTL/Lv_uORF/gene_pos_lv.txt", sep="")
positions = read.csv(gene_location_file_name, sep="\t", stringsAsFactors=F)
snp.pos = as.data.frame(fread(snps_location_file_name))
covariates=NULL
genes = genes_filtered_lv 

save(pheno, covariates, positions, snp.pos, genes, genotypes_alleles, genotypes_alleles, phenotype.file, genotypes.file, covariates.file,snps_location_file_name, gene_location_file_name, file = paste(res_dir, "Input_MatrixEQTL/Lv_uORF/Lv_Ribo_uORF_input.RData", sep=""))


load(paste(res_dir, "Input_MatrixEQTL/Lv_uORF/Lv_Ribo_uORF_input.RData", sep=""))

actual.eqtls_uORF_lv = eqtl_new(phenotype.file, genotypes.file, covariates.file, positions, snp.pos, paste(res_dir, "Results_MatrixEQTL/20180228_Lv_uORF/Lv_Ribo_uORF", sep=""), compute.all=T, only_cis=T)
save(actual.eqtls_uORF_lv, file=paste(res_dir, "Results_MatrixEQTL/20180228_Lv_uORF/eqtls.RData", sep=""))

#Liver

pheno = uORF_ribo_li
phenotype.file = paste(res_dir, "Input_MatrixEQTL/Liver_uORF/phenotypes.txt", sep="")

qtl_dir = "/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/"

genotypes.file =  paste(qtl_dir, "Input/sdp_genotypes.txt", sep="")
covariates.file = NULL
snps_location_file_name = paste(qtl_dir, "Input/snp_pos.txt", sep="")
gene_location_file_name = paste(res_dir, "Input_MatrixEQTL/Liver_uORF/gene_pos_li.txt", sep="")
positions = read.csv(gene_location_file_name, sep="\t", stringsAsFactors=F)
snp.pos = as.data.frame(fread(snps_location_file_name))
covariates=NULL
genes = genes_filtered_li 

save(pheno, covariates, positions, snp.pos, genes, genotypes_alleles, genotypes_alleles, phenotype.file, genotypes.file, covariates.file,snps_location_file_name, gene_location_file_name, file = paste(res_dir, "Input_MatrixEQTL/Liver_uORF/Liver_Ribo_uORF_input.RData", sep=""))


load(paste(res_dir, "Input_MatrixEQTL/Liver_uORF/Liver_Ribo_uORF_input.RData", sep=""))

actual.eqtls_uORF_li = eqtl_new(phenotype.file, genotypes.file, covariates.file, positions, snp.pos, paste(res_dir, "Results_MatrixEQTL/20180228_Liver_uORF/Liver_Ribo_uORF", sep=""), compute.all=T, only_cis=T)
save(actual.eqtls_uORF_li, file=paste(res_dir, "Results_MatrixEQTL/20180228_Liver_uORF/eqtls.RData", sep=""))

uORF_lv_sign = actual.eqtls_uORF_lv[which(actual.eqtls_uORF_lv$FDR < 0.05),] # 46 uORFs in 27 genes
uORF_li_sign = actual.eqtls_uORF_li[which(actual.eqtls_uORF_li$FDR < 0.05),] # 18 uORFs in 13 genes



load("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/permutations/All_tables_permutations.RData")
polyA_lv_with_uORF_QTL = data.frame(lv_polyA_local[which(lv_polyA_local$SNP %in% uORF_lv_sign$SNP),])
polyA_lv_with_uORF_QTL = data.frame(polyA_lv_with_uORF_QTL, uORF_lv_sign[match(polyA_lv_with_uORF_QTL$SNP, uORF_lv_sign$SNP),])
polyA_lv_with_uORF_QTL = polyA_lv_with_uORF_QTL[which(polyA_lv_with_uORF_QTL$gene_name == polyA_lv_with_uORF_QTL$gene_name.1),] #14 --> to check

ribo_lv_with_uORF_QTL = data.frame(lv_ribo_local[which(lv_ribo_local$SNP %in% uORF_lv_sign$SNP),])
ribo_lv_with_uORF_QTL = data.frame(ribo_lv_with_uORF_QTL, uORF_lv_sign[match(ribo_lv_with_uORF_QTL$SNP, uORF_lv_sign$SNP),])
ribo_lv_with_uORF_QTL = ribo_lv_with_uORF_QTL[which(ribo_lv_with_uORF_QTL$gene_name == ribo_lv_with_uORF_QTL$gene_name.1),]
#5 --> to check

te_lv_with_uORF_QTL = data.frame(lv_te_local[which(lv_te_local$SNP %in% uORF_lv_sign$SNP),])
te_lv_with_uORF_QTL = data.frame(te_lv_with_uORF_QTL, uORF_lv_sign[match(te_lv_with_uORF_QTL$SNP, uORF_lv_sign$SNP),])
te_lv_with_uORF_QTL = te_lv_with_uORF_QTL[which(te_lv_with_uORF_QTL$gene_name == te_lv_with_uORF_QTL$gene_name.1),] #2



polyA_li_with_uORF_QTL = data.frame(li_polyA_local[which(li_polyA_local$SNP %in% uORF_li_sign$SNP),])
polyA_li_with_uORF_QTL = data.frame(polyA_li_with_uORF_QTL, uORF_li_sign[match(polyA_li_with_uORF_QTL$SNP, uORF_li_sign$SNP),])
polyA_li_with_uORF_QTL = polyA_li_with_uORF_QTL[which(polyA_li_with_uORF_QTL$gene_name == polyA_li_with_uORF_QTL$gene_name.1),] #4 --> to check

ribo_li_with_uORF_QTL = data.frame(li_ribo_local[which(li_ribo_local$SNP %in% uORF_li_sign$SNP),])
ribo_li_with_uORF_QTL = data.frame(ribo_li_with_uORF_QTL, uORF_li_sign[match(ribo_li_with_uORF_QTL$SNP, uORF_li_sign$SNP),])
ribo_li_with_uORF_QTL = ribo_li_with_uORF_QTL[which(ribo_li_with_uORF_QTL$gene_name == ribo_li_with_uORF_QTL$gene_name.1),]
#5 --> to check

te_li_with_uORF_QTL = data.frame(li_te_local[which(li_te_local$SNP %in% uORF_li_sign$SNP),])
te_li_with_uORF_QTL = data.frame(te_li_with_uORF_QTL, uORF_li_sign[match(te_li_with_uORF_QTL$SNP, uORF_li_sign$SNP),])
te_li_with_uORF_QTL = te_li_with_uORF_QTL[which(te_li_with_uORF_QTL$gene_name == te_li_with_uORF_QTL$gene_name.1),] #n=0

save(polyA_lv_with_uORF_QTL, ribo_lv_with_uORF_QTL, te_lv_with_uORF_QTL, polyA_li_with_uORF_QTL, ribo_li_with_uORF_QTL, te_li_with_uORF_QTL, file = "/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment/significant_uORF_QTLs_overlapQTLs.RData")

library(ggplot2)


load(file = "/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment/significant_uORF_QTLs_overlapQTLs.RData")



res_dir="/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment/"
load(file=paste(res_dir, "Results_MatrixEQTL/20180228_Liver_uORF/eqtls.RData", sep=""))
load(file=paste(res_dir, "Results_MatrixEQTL/20180228_Lv_uORF/eqtls.RData", sep=""))
write.table(actual.eqtls_uORF_lv, file = "Desktop/Manuscript/uORF_QTL_results_lv.txt", sep="\t", col.names=T, row.names = F, quote = F)
write.table(actual.eqtls_uORF_li, file = "Desktop/Manuscript/uORF_QTL_results_liver.txt", sep="\t", col.names=T, row.names = F, quote = F)

