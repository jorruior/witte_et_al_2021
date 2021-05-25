##sort all qtl results by gene and sdp position
library(GenomicRanges)
library(eQTLpipeline)
source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/SeqQTL_oldCluster.R")



dir="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Result_tables/"
genes_names = gtf2GR("/data/huebner2/DATASETS/GENOMES/Ensembl_Genomes/Ens82/rn6/annotation/Rattus_norvegicus.Rnor_6.0.82.gtf", c("gene_id", "gene_name", "gene_biotype"))
genes_names = genes_names[values(genes_names)[,"type"] == "gene"]
genes_names = genes_names[which(seqnames(genes_names) %in% c(1:20, "X")),]

genotype_info = read.table("~/work/new_data_RI_translational_efficiency/Franzi/results/20160726/Input_MatrixEQTL/snp_pos.txt", sep="\t", header=T, stringsAsFactors=F) 
genotypes = read.table("~/work/RI_translational_efficiency/data/20160318/RN6_SDPs/HXBBXH-SDP_RN60_LDB-genome.txt",  sep="\t", comment.char = "#", header=F, stringsAsFactors=F, colClasses="character", check.names=F)
colnames(genotypes) = c("ID", "RI_geno", "SDP_local")
genotypes = genotypes[match(genotype_info$snp_id, genotypes$ID),]

local_sdps = read.table("/home/fkreuch/work/RI_translational_efficiency/data/20160318/RN6_SDPs/HXBBXH-SDP_RN60_LDB-local.txt", sep="\t", comment.char = "#", header=F, stringsAsFactors=F, colClasses="character", check.names=F)
colnames(local_sdps) = c("local_ID", "chr", "start", "end", "Marker_RN60_CSV", "Marker_count")
rownames(local_sdps) = local_sdps$local_ID
local_sdps = local_sdps[,c(1:4)]
SDP_index = apply(local_sdps, 1, function(x) grep(x[1], genotypes$SDP_local))
SDP_index = SDP_index[sapply(SDP_index, function(x) length(x)!=0)]
local_sdps = data.frame(local_sdps[match(names(SDP_index), local_sdps$local_ID),], genotypes[as.numeric(SDP_index), ])


local_SDPs_ranges = GRanges(seqnames = local_sdps$chr, ranges=IRanges(start=as.numeric(local_sdps$start), end = as.numeric(local_sdps$end)), global_ID = local_sdps$ID, local_ID = local_sdps$local_ID, add_local_IDs = local_sdps$SDP_local)
overlap_sdp_gene = findOverlaps(genes_names, local_SDPs_ranges, type="within")
genes_names$SDP_local = "noSDP"
genes_names[queryHits(overlap_sdp_gene)]$SDP_local = local_SDPs_ranges[subjectHits(overlap_sdp_gene)]$local_ID
genes_names$addSDP_local = "noSDP"
genes_names[queryHits(overlap_sdp_gene)]$addSDP_local = local_SDPs_ranges[subjectHits(overlap_sdp_gene)]$add_local_IDs

genes_names$SDP_global = "noSDP"
genes_names[queryHits(overlap_sdp_gene)]$SDP_global = local_SDPs_ranges[subjectHits(overlap_sdp_gene)]$global_ID



#lv polyA
load(file=paste(dir, "20180213_Lv_polyA/eqtls.RData", sep=""))
eQTL_results_lv_polyA = as.data.frame(rbind(actual.eqtls[["cis"]], actual.eqtls[["trans"]]))
eQTL_results_lv_polyA$gene_contSDP = genes_names[match(eQTL_results_lv_polyA$gene, genes_names$gene_id), ]$SDP_global
eQTL_results_lv_polyA$SDP_local = genotypes[match(eQTL_results_lv_polyA$SNP, genotypes$ID), "SDP_local"]
eQTL_results_lv_polyA[which(eQTL_results_lv_polyA[,9] == "X"), 9] = 21
eQTL_results_lv_polyA = eQTL_results_lv_polyA[order(eQTL_results_lv_polyA$gene, eQTL_results_lv_polyA$snp_pos),]
save(eQTL_results_lv_polyA, file=paste(dir, "20180213_Lv_polyA/eqtls_sorted.RData", sep=""))

#liver polyA
load(file=paste(dir, "20180213_Liver_polyA/eqtls.RData", sep=""))
eQTL_results_liver_polyA = as.data.frame(rbind(actual.eqtls[["cis"]], actual.eqtls[["trans"]]))
eQTL_results_liver_polyA$gene_contSDP = genes_names[match(eQTL_results_liver_polyA$gene, genes_names$gene_id), ]$SDP_global
eQTL_results_liver_polyA$SDP_local = genotypes[match(eQTL_results_liver_polyA$SNP, genotypes$ID), "SDP_local"]
eQTL_results_liver_polyA[which(eQTL_results_liver_polyA[,9] == "X"), 9] = 21
eQTL_results_liver_polyA = eQTL_results_liver_polyA[order(eQTL_results_liver_polyA$gene, eQTL_results_liver_polyA$snp_pos),]
save(eQTL_results_liver_polyA, file=paste(dir, "20180213_Liver_polyA/eqtls_sorted.RData", sep=""))

#lv ribo
load(file=paste(dir, "20180213_Lv_Ribo/eqtls.RData", sep=""))
eQTL_results_lv_ribo = as.data.frame(rbind(actual.eqtls[["cis"]], actual.eqtls[["trans"]]))
eQTL_results_lv_ribo$gene_contSDP = genes_names[match(eQTL_results_lv_ribo$gene, genes_names$gene_id), ]$SDP_global
eQTL_results_lv_ribo$SDP_local = genotypes[match(eQTL_results_lv_ribo$SNP, genotypes$ID), "SDP_local"]
eQTL_results_lv_ribo[which(eQTL_results_lv_ribo[,9] == "X"), 9] = 21
eQTL_results_lv_ribo = eQTL_results_lv_ribo[order(eQTL_results_lv_ribo$gene, eQTL_results_lv_ribo$snp_pos),]
save(eQTL_results_lv_ribo, file=paste(dir, "20180213_Lv_Ribo/eqtls_sorted.RData", sep=""))

#liver ribo
load(file=paste(dir, "20180213_Liver_Ribo/eqtls.RData", sep=""))
eQTL_results_liver_ribo = as.data.frame(rbind(actual.eqtls[["cis"]], actual.eqtls[["trans"]]))
eQTL_results_liver_ribo$gene_contSDP = genes_names[match(eQTL_results_liver_ribo$gene, genes_names$gene_id), ]$SDP_global
eQTL_results_liver_ribo$SDP_local = genotypes[match(eQTL_results_liver_ribo$SNP, genotypes$ID), "SDP_local"]
eQTL_results_liver_ribo[which(eQTL_results_liver_ribo[,9] == "X"), 9] = 21
eQTL_results_liver_ribo = eQTL_results_liver_ribo[order(eQTL_results_liver_ribo$gene, eQTL_results_liver_ribo$snp_pos),]
save(eQTL_results_liver_ribo, file=paste(dir, "20180213_Liver_Ribo/eqtls_sorted.RData", sep=""))

#lv TE
load(file=paste(dir, "20180213_Lv_Ribo_RNA_resid/eqtls.RData", sep=""))
eQTL_results_lv_riboRnaResid = as.data.frame(rbind(actual.eqtls[["cis"]], actual.eqtls[["trans"]]))
eQTL_results_lv_riboRnaResid$gene_contSDP = genes_names[match(eQTL_results_lv_riboRnaResid$gene, genes_names$gene_id), ]$SDP_global
eQTL_results_lv_riboRnaResid$SDP_local = genotypes[match(eQTL_results_lv_riboRnaResid$SNP, genotypes$ID), "SDP_local"]
eQTL_results_lv_riboRnaResid[which(eQTL_results_lv_riboRnaResid[,9] == "X"), 9] = 21
eQTL_results_lv_riboRnaResid = eQTL_results_lv_riboRnaResid[order(eQTL_results_lv_riboRnaResid$gene, eQTL_results_lv_riboRnaResid$snp_pos),]
save(eQTL_results_lv_riboRnaResid, file=paste(dir, "20180213_Lv_Ribo_RNA_resid/eqtls_sorted.RData", sep=""))

#liver TE
load(file=paste(dir, "20180403_Liver_Ribo_RNA_resid/eqtls.RData", sep=""))
eQTL_results_liver_riboRnaResid = as.data.frame(rbind(actual.eqtls[["cis"]], actual.eqtls[["trans"]]))
eQTL_results_liver_riboRnaResid$gene_contSDP = genes_names[match(eQTL_results_liver_riboRnaResid$gene, genes_names$gene_id), ]$SDP_global
eQTL_results_liver_riboRnaResid$SDP_local = genotypes[match(eQTL_results_liver_riboRnaResid$SNP, genotypes$ID), "SDP_local"]
eQTL_results_liver_riboRnaResid[which(eQTL_results_liver_riboRnaResid[,9] == "X"), 9] = 21
eQTL_results_liver_riboRnaResid = eQTL_results_liver_riboRnaResid[order(eQTL_results_liver_riboRnaResid$gene, eQTL_results_liver_riboRnaResid$snp_pos),]
save(eQTL_results_liver_riboRnaResid, file=paste(dir, "20180403_Liver_Ribo_RNA_resid/eqtls_sorted.RData", sep=""))





