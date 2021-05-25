#library(eQTLpipeline)
library("MatrixEQTL", lib.loc="/home/jruizor/lib/")
library(data.table)
library(GenomicRanges)
library(grid)

source("/data/huebner2/ANALYSES/20151127_qtl_RI_rn5/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/eqtl_trans_kinship.R")
source("/data/huebner2/ANALYSES/20151127_qtl_RI_rn5/RI_translational_efficiency/scripts/R//eQTLpipeline_functions/SeqQTL_oldCluster.R")

#sva = read.table("Desktop/sva_all.txt", sep="\t", header=T)

res_dir="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/"
data_dir="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/"
load(file=paste(data_dir, "preprocessed_data.RData", sep=""))

sample_names = c("BXH02", "BXH03", "BXH05", "HXB01", "HXB02", "HXB03", "HXB04", "HXB05", "HXB07", "HXB10", "HXB13", "HXB15", "HXB17", "HXB18", "HXB20", "HXB21", "HXB22", "HXB23", "HXB24", "HXB25", "HXB27", "HXB29", "HXB31", "BXH06", "BXH08", "BXH09", "BXH10", "BXH11", "BXH12", "BXH13") 

#Phenotype Input file for eQTL Mapping
log2_polyA_lv = log2_polyA_lv[,match(sample_names, colnames(log2_polyA_lv))]
phenotype.file1 = paste(res_dir, "Input/Lv_polyA/phenotypes.txt", sep="")
#write.table(log2_polyA_lv, phenotype.file1, sep="\t", quote=F, row.names=T, col.names=T)

log2_polyA_li = log2_polyA_li[,match(sample_names, colnames(log2_polyA_li))]
phenotype.file2 = paste(res_dir, "Input/Liver_polyA/phenotypes.txt", sep="")
#write.table(log2_polyA_li, phenotype.file2, sep="\t", quote=F, row.names=T, col.names=T)

log2_ribo_lv = log2_ribo_lv[,match(sample_names, colnames(log2_ribo_lv))]
phenotype.file3 = paste(res_dir, "Input/Lv_Ribo/phenotypes.txt", sep="")
#write.table(log2_ribo_lv, phenotype.file3, sep="\t", quote=F, row.names=T, col.names=T)

colnames(log2_ribo_li) = gsub("_4", "", colnames(log2_ribo_li))
log2_ribo_li = log2_ribo_li[,match(sample_names, colnames(log2_ribo_li))]
phenotype.file4 = paste(res_dir, "Input/Liver_Ribo/phenotypes.txt", sep="")
#write.table(log2_ribo_li, phenotype.file4, sep="\t", quote=F, row.names=T, col.names=T)

print("Data loaded")

#Make a gene position file of the expressed genes that need to be tested
genes_filtered_lv = data.frame(gene_id = genes_names[match(rownames(log2_ribo_lv), genes_names$gene_id), ]$gene_id, chrom = seqnames(genes_names[match(rownames(log2_ribo_lv), genes_names$gene_id), ]), start = start(genes_names[match(rownames(log2_ribo_lv), genes_names$gene_id), ]), end = end(genes_names[match(rownames(log2_ribo_lv), genes_names$gene_id), ]), gene_name=genes_names[match(rownames(log2_ribo_lv), genes_names$gene_id), ]$gene_name, gene_type=genes_names[match(rownames(log2_ribo_lv), genes_names$gene_id), ]$type)

#write.table(genes_filtered_lv, file = paste(res_dir, "Input/gene_pos_lv.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)

genes_filtered_li = data.frame(gene_id = genes_names[match(rownames(log2_ribo_li), genes_names$gene_id), ]$gene_id, chrom = seqnames(genes_names[match(rownames(log2_ribo_li), genes_names$gene_id), ]), start = start(genes_names[match(rownames(log2_ribo_li), genes_names$gene_id), ]), end = end(genes_names[match(rownames(log2_ribo_li), genes_names$gene_id), ]), gene_name=genes_names[match(rownames(log2_ribo_li), genes_names$gene_id), ]$gene_name, gene_type=genes_names[match(rownames(log2_ribo_li), genes_names$gene_id), ]$type)
#write.table(genes_filtered_li, file = paste(res_dir, "Input/gene_pos_liver.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)


#load genotypes to convert Olis file into a matrix EQTL suitable format
genotypes = read.table("/data/huebner2/ANALYSES/20151127_qtl_RI_rn5/RI_translational_efficiency/data/20160318/RN6_SDPs/HXBBXH-SDP_RN60_LDB-genome.txt",  sep="\t", comment.char = "#", header=F, stringsAsFactors=F, colClasses="character", check.names=F)
colnames(genotypes) = c("ID", "RI_geno", "SDP_local")

####kick out single-marker SDPs
local_sdps = read.table("/data/huebner2/ANALYSES/20151127_qtl_RI_rn5/RI_translational_efficiency/data/20160318/RN6_SDPs/HXBBXH-SDP_RN60_LDB-local.txt", sep="\t", comment.char = "#", header=F, stringsAsFactors=F, colClasses="character", check.names=F)
colnames(local_sdps) = c("local_ID", "chr", "start", "end", "Marker_RN60_CSV", "Marker_count")

identify_single_marker = which(local_sdps$Marker_count == 1)

sdp_before = local_sdps[identify_single_marker-1,"local_ID"]
sdp_after = local_sdps[identify_single_marker+1,"local_ID"]
sdp_after = sdp_after[-1]
id_before = sapply(sdp_before, function(x) genotypes[grep(x, genotypes$SDP_local),"ID"]) 
id_after = sapply(sdp_after, function(x) genotypes[grep(x, genotypes$SDP_local),"ID"]) 
sdps_to_delete = id_before[which(id_before == id_after)]

##only 1685 SDPs left
genotypes = genotypes[-which(genotypes$ID %in% sdps_to_delete),]
#"snp_pos" in Matthias description
sdp_info = data.frame(snp_id = genotypes$ID, chrom = sapply(genotypes$ID, function(x){as.numeric(substr(x, 6,7))}), snp_pos = sapply(genotypes$ID, function(x){as.numeric(substr(x, 8,nchar(x)))}))
#write.table(sdp_info,file = paste(res_dir, "Input/snp_pos.txt", sep=""), sep="\t", quote=F, row.names=F)

#dosage file in Mathias description
genotypes_alleles = data.frame(apply(genotypes, 1, function(x){strsplit(x[2], split="")}))
genotypes_alleles = t(genotypes_alleles[c(7:36),])
colnames(genotypes_alleles) = sample_names
rownames(genotypes_alleles) = sdp_info$snp_id
genotypes_alleles = genotypes_alleles[,match(sample_names, colnames(genotypes_alleles))]
#write.table(genotypes_alleles, file = paste(res_dir, "Input/sdp_genotypes.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=T)

#covar_lv_rna = read.table(paste(res_dir, "Input/covariates_peer_eQTL_lv.txt", sep=""), sep=" ", stringsAsFactors = F, header = T)
#covar_lv_rna = t(data.frame(id = colnames(covar_lv_rna), t(covar_lv_rna)))
#write.table(covar_lv_rna, file = paste(res_dir, "Input/covariates_sva_eQTL_lv_mod.txt", sep=""), sep="\t", col.names = F, row.names=T, quote = F)
#covar_lv_ribo = read.table(paste(res_dir, "Input/covariates_peer_rQTL_lv.txt", sep=""), sep=" ", stringsAsFactors = F, header = T)
#covar_lv_ribo = t(data.frame(id = colnames(covar_lv_ribo), t(covar_lv_ribo)))
#write.table(covar_lv_ribo, file = paste(res_dir, "Input/covariates_sva_rQTL_lv_mod.txt", sep=""), sep="\t", col.names = F, row.names=T, quote = F)
#covar_liver_rna = read.table(paste(res_dir, "Input/covariates_peer_eQTL_liver.txt", sep=""), sep=" ", stringsAsFactors = F, header = T)
#covar_liver_rna = t(data.frame(id = colnames(covar_liver_rna), t(covar_liver_rna)))
#write.table(covar_liver_rna, file = paste(res_dir, "Input/covariates_sva_eQTL_liver_mod.txt", sep=""), sep="\t", col.names = F, row.names=T, quote = F)
#covar_liver_ribo = read.table(paste(res_dir, "Input/covariates_peer_rQTL_liver.txt", sep=""), sep=" ", stringsAsFactors = F, header = T)
#covar_liver_ribo = t(data.frame(id = colnames(covar_liver_ribo), t(covar_liver_ribo)))
#write.table(covar_liver_ribo, file = paste(res_dir, "Input/covariates_peer_rQTL_liver_mod.txt", sep=""), sep="\t", col.names = F, row.names=T, quote = F)

print("Start mapping")

################## eQTL-Mapping #################
print("polyA lv")

### define input
phenotype.file = paste(res_dir, "Input/Lv_polyA/phenotypes.txt", sep="")
genotypes.file =  paste(res_dir, "Input/sdp_genotypes.txt", sep="")
covariates.file =  paste(res_dir, "Input/eQTL_lv_peer.tsv", sep="")
snps_location_file_name = paste(res_dir, "Input/snp_pos.txt", sep="")
gene_location_file_name = paste(res_dir, "Input/gene_pos_lv.txt", sep="")
positions = read.csv(gene_location_file_name, sep="\t", stringsAsFactors=F)
snp.pos = as.data.frame(fread(snps_location_file_name))
pheno=log2_polyA_lv
#covariates=covar_sva
genes = genes_filtered_lv 
save(pheno, positions, snp.pos, genes, genotypes_alleles, genotypes_alleles, phenotype.file, genotypes.file, covariates.file,snps_location_file_name, gene_location_file_name, file = paste(res_dir, "Input/Lv_polyA/Lv_polyA_input_peer.RData", sep=""))

print("polyA liver")

### define input
phenotype.file = paste(res_dir, "Input/Liver_polyA/phenotypes.txt", sep="")
genotypes.file =  paste(res_dir, "Input/sdp_genotypes.txt", sep="")
covariates.file =  paste(res_dir, "Input/eQTL_liv_peer.tsv", sep="")
snps_location_file_name = paste(res_dir, "Input/snp_pos.txt", sep="")
gene_location_file_name = paste(res_dir, "Input/gene_pos_liver.txt", sep="")
positions = read.csv(gene_location_file_name, sep="\t", stringsAsFactors=F)
snp.pos = as.data.frame(fread(snps_location_file_name))
pheno=log2_polyA_li
#covariates=covar_sva
genes = genes_filtered_li 
save(pheno, positions, snp.pos, genes, genotypes_alleles, phenotype.file, genotypes.file, covariates.file,snps_location_file_name, gene_location_file_name, file = paste(res_dir, "Input/Liver_polyA/Liver_polyA_input_peer.RData", sep=""))

print("ribo lv")

### define input
phenotype.file = paste(res_dir, "Input/Lv_Ribo/phenotypes.txt", sep="")
genotypes.file =  paste(res_dir, "Input/sdp_genotypes.txt", sep="")
covariates.file =  paste(res_dir, "Input/rQTL_lv_peer.tsv", sep="")
snps_location_file_name = paste(res_dir, "Input/snp_pos.txt", sep="")
gene_location_file_name = paste(res_dir, "Input/gene_pos_lv.txt", sep="")
positions = read.csv(gene_location_file_name, sep="\t", stringsAsFactors=F)
snp.pos = as.data.frame(fread(snps_location_file_name))
pheno=log2_ribo_lv
#covariates=covar_sva
genes = genes_filtered_lv 
save(pheno, positions, snp.pos, genes, genotypes_alleles, phenotype.file, genotypes.file, covariates.file,snps_location_file_name, gene_location_file_name, file = paste(res_dir, "Input/Lv_Ribo/Lv_ribo_input_peer.RData", sep=""))

print("ribo liver")

### define input
phenotype.file = paste(res_dir, "Input/Liver_Ribo/phenotypes.txt", sep="")
genotypes.file =  paste(res_dir, "Input/sdp_genotypes.txt", sep="")
covariates.file1 =  paste(res_dir, "Input/eQTL_lv_peer.tsv.5cov.tsv", sep="")
covariates.file2 =  paste(res_dir, "Input/eQTL_liv_peer.tsv", sep="")
covariates.file3 =  paste(res_dir, "Input/rQTL_lv_peer.tsv", sep="")
covariates.file4 =  paste(res_dir, "Input/rQTL_liv_peer.tsv", sep="")
covariates.file5 =  paste(res_dir, "Input/teQTL_lv_peer.tsv", sep="")
covariates.file6 =  paste(res_dir, "Input/teQTL_liv_peer.tsv", sep="")

snps_location_file_name = paste(res_dir, "Input/snp_pos.txt", sep="")
gene_location_file_name = paste(res_dir, "Input/gene_pos_liver.txt", sep="")
positions = read.csv(gene_location_file_name, sep="\t", stringsAsFactors=F)
snp.pos = as.data.frame(fread(snps_location_file_name))
pheno=log2_ribo_li
#covariates=covar_sva
genes = genes_filtered_li 
save(pheno, positions, snp.pos, genes, genotypes_alleles, phenotype.file, genotypes.file, covariates.file,snps_location_file_name, gene_location_file_name, file = paste(res_dir, "Input/Liver_Ribo/Liver_ribo_input_peer.RData", sep=""))


#library(openxlsx)
print("compute eQTLs")
load(paste(res_dir, "Input/Lv_polyA/Lv_polyA_input_peer.RData", sep=""))
actual.eqtls = eqtl_new(phenotype.file, genotypes.file, covariates.file1, positions, snp.pos, paste(res_dir, "Result_tables/20201117_Lv_polyA_peer/Lv_polyA", sep=""), compute.all=T, only_cis=F)
save(actual.eqtls, file = paste(res_dir, "Result_tables/20201117_Lv_polyA_peer/eqtls.RData", sep=""))

##EXTRA
print("EXTRA")
phenotype.file = paste(res_dir, "Input/Lv_polyA/phenotypes_sig.txt", sep="")
genotypes.file =  paste(res_dir, "Input/Lv_polyA/sdp_genotypes_sig.txt", sep="")
#positions = read.table(paste(res_dir, "Input/Lv_polyA/gene_pos_sig.txt", sep=""))
#snp.pos = as.data.frame(fread(paste(res_dir, "Input/Lv_polyA/snp_pos_sig.txt", sep="")))
actual.eqtls = eqtl_new(phenotype.file, genotypes.file, covariates.file1, positions, snp.pos, paste(res_dir, "Result_tables/20201117_Lv_polyA_sig_peer/Lv_polyA_sig", sep=""), compute.all=T, only_cis=F)
save(actual.eqtls, file = paste(res_dir, "Result_tables/20201117_Lv_polyA_sig_peer/eqtls_sig.RData", sep=""))
##

load(paste(res_dir, "Input/Liver_polyA/Liver_polyA_input_peer.RData", sep=""))
actual.eqtls = eqtl_new(phenotype.file, genotypes.file, covariates.file1, positions, snp.pos, paste(res_dir, "Result_tables/20201119_Liver_polyA_peer/Liver_polyA", sep=""), compute.all=T, only_cis=F)
save(actual.eqtls, file=paste(res_dir, "Result_tables/20201119_Liver_polyA_peer/eqtls.RData", sep=""))

load(paste(res_dir, "Input/Lv_Ribo/Lv_ribo_input_peer.RData", sep=""))
actual.eqtls = eqtl_new(phenotype.file, genotypes.file, covariates.file1, positions, snp.pos, paste(res_dir, "Result_tables/20201117_Lv_Ribo_peer/Lv_ribo", sep=""), compute.all=T, only_cis=F)
save(actual.eqtls, file=paste(res_dir, "Result_tables/20201117_Lv_Ribo_peer/eqtls.RData", sep=""))

load(paste(res_dir, "Input/Liver_Ribo/Liver_ribo_input.RData", sep=""))
actual.eqtls = eqtl_new(phenotype.file, genotypes.file, covariates.file1, positions, snp.pos, paste(res_dir, "Result_tables/20201119_Liver_Ribo_peer/Liver_ribo", sep=""), compute.all=T, only_cis=F)
save(actual.eqtls, file=paste(res_dir, "Result_tables/20201119_Liver_Ribo_peer/eqtls.RData", sep=""))

########### TE - use residuals of lm(ribo~rna) ########### 


#covar_lv_te = read.table(paste(res_dir, "Input/covar", sep=""), sep=" ", stringsAsFactors = F, header = T)
#covar_lv_te = t(data.frame(id = colnames(covar_lv_te), t(covar_lv_te)))
#write.table(covar_lv_te, file = paste(res_dir, "Input/covariates_peer_teQTL_lv_mod.txt", sep=""), sep="\t", col.names = F, row.names=T, quote = F)


pheno=resid(lm(as.matrix(log2_ribo_lv)~as.matrix(log2_polyA_lv)))
pheno1=log2_ribo_lv-log2_polyA_lv
TE=pheno

print("TEqtl lv")

phenotype.file = paste(res_dir, "Input/Lv_Ribo_RNA_resid/phenotypes.txt", sep="")
write.table(pheno, phenotype.file, sep="\t", quote=F, row.names=T, col.names=T)
genotypes.file =  paste(res_dir, "Input/sdp_genotypes.txt", sep="")
covariates.file =  paste(res_dir, "Input/eQTL_lv_peer.tsv", sep="")
snps_location_file_name = paste(res_dir, "Input/snp_pos.txt", sep="")
gene_location_file_name = paste(res_dir, "Input/gene_pos_lv.txt", sep="")
positions = read.csv(gene_location_file_name, sep="\t", stringsAsFactors=F)
snp.pos = as.data.frame(fread(snps_location_file_name))
#covariates=covar_sva
genes = genes_filtered_lv 
save(pheno, positions, snp.pos, genes, genotypes_alleles, genotypes_alleles, phenotype.file, genotypes.file, covariates.file,snps_location_file_name, gene_location_file_name, file = paste(res_dir, "Input/Lv_Ribo_RNA_resid/Lv_Ribo_RNA_resid_input_peer.RData", sep=""))


load(paste(res_dir, "Input/Lv_Ribo_RNA_resid/Lv_Ribo_RNA_resid_input_peer.RData", sep=""))
actual.eqtls = eqtl_new(phenotype.file, genotypes.file, covariates.file, positions, snp.pos, paste(res_dir, "Result_tables/20201117_Lv_Ribo_RNA_resid_peer/Lv_Ribo_RNA_resid_", sep=""), compute.all=T, only_cis=F)
save(actual.eqtls, file=paste(res_dir, "Result_tables/20201117_Lv_Ribo_RNA_resid_peer/eqtls.RData", sep=""))


print("TEqtl liver")

### define input liver
pheno=resid(lm(as.matrix(log2_ribo_li)~as.matrix(log2_polyA_li)))
phenotype.file = paste(res_dir, "Input/Liver_Ribo_RNA_resid/phenotypes.txt", sep="")
write.table(pheno, phenotype.file, sep="\t", quote=F, row.names=T, col.names=T)

genotypes.file =  paste(res_dir, "Input/sdp_genotypes.txt", sep="")
covariates.file = paste(res_dir, "Input/eQTL_liv_peer.tsv", sep="")
snps_location_file_name = paste(res_dir, "Input/snp_pos.txt", sep="")
gene_location_file_name = paste(res_dir, "Input/gene_pos_liver.txt", sep="")
positions = read.csv(gene_location_file_name, sep="\t", stringsAsFactors=F)
snp.pos = as.data.frame(fread(snps_location_file_name))
#covariates=NULL
genes = genes_filtered_li 
save(pheno, positions, snp.pos, genes, genotypes_alleles, genotypes_alleles, phenotype.file, genotypes.file, covariates.file,snps_location_file_name, gene_location_file_name, file = paste(res_dir, "Input/Liver_Ribo_RNA_resid/Liver_Ribo_RNA_resid_input_peer.RData", sep=""))


load(paste(res_dir, "Input/Liver_Ribo_RNA_resid/Liver_Ribo_RNA_resid_input.RData", sep=""))
actual.eqtls = eqtl_new(phenotype.file, genotypes.file, covariates.file, positions, snp.pos, paste(res_dir, "Result_tables/20180403_Liver_Ribo_RNA_resid_peer/Liver_Ribo_RNA_resid_", sep=""), compute.all=T, only_cis=F)
save(actual.eqtls, file=paste(res_dir, "Result_tables/20180403_Liver_Ribo_RNA_resid_peer/eqtls.RData", sep=""))



