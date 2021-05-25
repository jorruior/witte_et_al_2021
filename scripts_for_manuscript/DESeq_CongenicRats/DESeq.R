library(eQTLpipeline)
library(preprocessCore)
library(GenomicRanges)
library(DESeq2)
source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/SeqQTL_oldCluster.R")
source("~/work/HS_heart/R/functions/vif_function.R")
source("~/work/HS_heart/R/functions/rahmann_scaling.R")


load(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/filtered_genesets_ribotaper_fpkm.RData")

res_dir="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/DESeq_congenics/"

genes_names = gtf2GR("/data/huebner2/DATASETS/GENOMES/Ensembl_Genomes/Ens82/rn6/annotation/Rattus_norvegicus.Rnor_6.0.82.gtf", c("gene_id", "gene_name", "gene_biotype"))
genes_names = genes_names[values(genes_names)[,"type"] == "gene"]
genes_names = genes_names[which(seqnames(genes_names) %in% c(1:20, "X")),]

ribo_raw_ref_lv = read.table(file = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/raw/ribo_cong_lv_proteincodingANDlincRNA_countMatrix_raw.txt", sep="\t", header=T)

ribo_raw_ref_lv = ribo_raw_ref_lv[which(rownames(ribo_raw_ref_lv) %in% genesset_lv), ]
colnames(ribo_raw_ref_lv) = gsub("_counts_in_CDS_stranded.txt", "", colnames(ribo_raw_ref_lv))
polyA_raw_ref_lv = read.table(file = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/raw/polyA_cong_lv_proteincodingANDlincRNA_countMatrix_raw.txt", sep="\t", header=T)

polyA_raw_ref_lv = polyA_raw_ref_lv[which(rownames(polyA_raw_ref_lv) %in% genesset_lv), ]
colnames(polyA_raw_ref_lv) = gsub("_counts_in_CDS_strandedreverse.txt", "", colnames(polyA_raw_ref_lv))

pool_sizeFactor <- estimateSizeFactorsForMatrix(cbind(polyA_raw_ref_lv,ribo_raw_ref_lv))
mrna_sizeFactor <- pool_sizeFactor[1:ncol(polyA_raw_ref_lv)]
rpf_sizeFactor <- pool_sizeFactor[(ncol(polyA_raw_ref_lv)+1):(ncol(polyA_raw_ref_lv)+ncol(ribo_raw_ref_lv))]


all_counts = cbind(ribo_raw_ref_lv)

colData_ribo = data.frame(strain=c(rep("BN", 5), rep("SHR", 5)), row.names = colnames(all_counts))


summary_ribo_lv = DESeqDataSetFromMatrix(countData=all_counts, colData_ribo ,design=~strain)
colData(summary_ribo_lv)$genome=factor(colData(summary_ribo_lv)$strain, levels=c("BN","SHR"))
dds_ribo_lv = summary_ribo_lv
sizeFactors(dds_ribo_lv) = rpf_sizeFactor
colnames(dds_ribo_lv) = colnames(all_counts)
dds_ribo_lv_res = DESeq(dds_ribo_lv)
diff_res_ribo_lv = results(dds_ribo_lv_res)

diff_res_ribo_lv = na.omit(diff_res_ribo_lv)
diff_res_ribo_lv = data.frame(diff_res_ribo_lv, genes_names[match(rownames(diff_res_ribo_lv), genes_names$gene_id),]$gene_name)

write.table(diff_res_ribo_lv[which(diff_res_ribo_lv$padj <= 0.05),], file = paste(res_dir, "ribo_deseq2_results_congenics_lv.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)

diff_res_ribo_lv$sign = (diff_res_ribo_lv$padj <= 0.05 & (diff_res_ribo_lv$log2FoldChange >= log2(1.2) | diff_res_ribo_lv$log2FoldChange <= log2(1/1.2)))



##trimmed PolyA


all_counts = cbind(polyA_raw_ref_lv)

colData_polyA = data.frame(strain=c(rep("BN", 5), rep("SHR", 5)), row.names = colnames(all_counts))

summary_polyA_lv = DESeqDataSetFromMatrix(countData=all_counts, colData_polyA ,design=~strain)
colData(summary_polyA_lv)$genome=factor(colData(summary_polyA_lv)$strain, levels=c("BN","SHR"))
dds_polyA_lv = summary_polyA_lv
sizeFactors(dds_polyA_lv) = mrna_sizeFactor

colnames(dds_polyA_lv) = colnames(all_counts)
dds_polyA_lv_res = DESeq(dds_polyA_lv)
diff_res_polyA_lv = results(dds_polyA_lv_res)

diff_res_polyA_lv = na.omit(diff_res_polyA_lv)
diff_res_polyA_lv = data.frame(diff_res_polyA_lv, genes_names[match(rownames(diff_res_polyA_lv), genes_names$gene_id),]$gene_name)

write.table(diff_res_polyA_lv[which(diff_res_polyA_lv$padj <= 0.05),], file = paste(res_dir, "polyA_deseq2_results_congenics_lv.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)

diff_res_polyA_lv$sign = (diff_res_polyA_lv$padj <= 0.05 & (diff_res_polyA_lv$log2FoldChange >= log2(1.2) | diff_res_polyA_lv$log2FoldChange <= log2(1/1.2)))


table(diff_res_polyA_lv$sign)
table(diff_res_ribo_lv$sign)


save(diff_res_polyA_lv, diff_res_ribo_lv, file = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/DESeq_congenics/result_deseq.RData")


#822 genes differentially translated and not on chr. 3
diff_res_ribo_lv$chr = seqnames(genes_names[match(rownames(diff_res_ribo_lv), genes_names$gene_id),])
ribo_notChr3 = diff_res_ribo_lv[which(diff_res_ribo_lv$chr !=3 & diff_res_ribo_lv$sign == T),]
diff_res_polyA_lv$chr = seqnames(genes_names[match(rownames(diff_res_polyA_lv), genes_names$gene_id),])
polyA_notChr3 = diff_res_polyA_lv[which(diff_res_polyA_lv$chr !=3 & diff_res_polyA_lv$sign == T),] #1433


#716 purely regulated on ribo level and not located on chr. 3
ribo_only_notChr3 = ribo_notChr3[setdiff(rownames(ribo_notChr3), rownames(polyA_notChr3)),]

library(gProfileR)
set_base_url("http://biit.cs.ut.ee/gprofiler_archive2/r1477_e82_eg29/web")
#716
ribo_only_notChr3_go = gprofiler(query =  unique(rownames(ribo_only_notChr3)), custom_bg = genesset_lv, organism = "rnorvegicus")

ribo_all_go = gprofiler(query =  unique(rownames(diff_res_ribo_lv[which(diff_res_ribo_lv$sign == T),])), custom_bg = genesset_lv, organism = "rnorvegicus")


polyA_all_go = gprofiler(query =  unique(rownames(diff_res_polyA_lv[which(diff_res_polyA_lv$sign == T),])), custom_bg = genesset_lv, organism = "rnorvegicus")

save(diff_res_polyA_lv, diff_res_ribo_lv, file = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/DESeq_congenics/result_deseq.RData")
save(ribo_only_notChr3_go, ribo_only_notChr3, file = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/DESeq_congenics/FC_plot_riboOnlyGenesNotLocatedChr3.RData")

load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/permutations/All_tables_permutations.RData")
matrixEQTL_distantRibo = gprofiler(query =  unique(lv_ribo_distant[which(lv_ribo_distant$SNP %in% c("SDPG_03006314843", "SDPG_03009805651", "SDPG_03003960268")),]$gene), custom_bg = genesset_lv, organism = "rnorvegicus")


###chr3 cluster ECM genes
ECM_congenics = unlist(strsplit(ribo_only_notChr3_go[83,"intersection"], ",")



#make a heatmap of all expressed genes (deseq normalized and scaled)

ribo_deseq_cong_lv = read.table(file = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/deseq_norm/polyA_cong_lv_Ribo_deseqNorm_RNAANDRIBO.txt", sep="\t", header=T)
polya_deseq_cong_lv = read.table(file = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/deseq_norm/polyA_cong_lv_RNA_deseqNorm_RNAANDRIBO.txt", sep="\t", header=T)


load(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/filtered_genesets_ribotaper_fpkm.RData")


hmcol = colorRampPalette(c("darkred", "white", "darkblue"))(n = 10)
library(gplots)

scaled_counts_ribo = t(scale(t(log2(ribo_deseq_cong_lv[genesset_lv,]+1))))
colnames(scaled_counts_ribo) = gsub("_Ri_counts_in_CDS_stranded.txt", "", colnames(scaled_counts_ribo))
scaled_counts_polyA = t(scale(t(log2(polya_deseq_cong_lv[genesset_lv,]+1))))
colnames(scaled_counts_polyA) = gsub("_mR_counts_in_CDS_strandedreverse.txt", "", colnames(scaled_counts_ribo))

means_scaled_counts_ribo = data.frame(endo = apply(scaled_counts_ribo[,1:5], 1, mean), ctrl = apply(scaled_counts_ribo[,6:10], 1, mean))
means_scaled_counts_polyA = data.frame(endo = apply(scaled_counts_polyA[,1:5], 1, mean), ctrl = apply(scaled_counts_polyA[,6:10], 1, mean))



pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/DESeq_congenics/scaled_counts.pdf", width = 15)
heatmap.2(scaled_counts_ribo, main = "Congenics Ribo" ,col=hmcol, trace="none", symbreaks=T,labRow=F,Colv="NA", dendrogram="none", margin=c(8,8))
heatmap.2(as.matrix(means_scaled_counts_ribo), main = "Congenics meanRibo" ,col=hmcol, trace="none", symbreaks=T,labRow=F,Colv="NA", dendrogram="none", margin=c(8,8))

heatmap.2(na.omit(scaled_counts_polyA), main = "Congenics PolyA" ,col=hmcol, trace="none", symbreaks=T,labRow=F,Colv="NA", dendrogram="none", margin=c(8,8))
heatmap.2(as.matrix(na.omit(means_scaled_counts_polyA)), main = "Congenics meanPolyA" ,col=hmcol, trace="none", symbreaks=T,labRow=F,Colv="NA", dendrogram="none", margin=c(8,8))

dev.off()
