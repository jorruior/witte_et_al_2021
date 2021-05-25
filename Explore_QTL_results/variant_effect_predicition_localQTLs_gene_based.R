
source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/SeqQTL_oldCluster.R")

library(GenomicRanges)
load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/permutations/All_tables_permutations.RData")

#variant effect prediction based on SNPs in local QTLs

all_locally_regulated_genes = unique(c(lv_polyA_local$gene, lv_ribo_local$gene,  lv_te_local$gene, li_polyA_local$gene, li_ribo_local$gene, li_te_local$gene))

genes_names = gtf2GR("/data/huebner2/DATASETS/GENOMES/Ensembl_Genomes/Ens82/rn6/annotation/Rattus_norvegicus.Rnor_6.0.82.gtf", c("gene_id", "gene_name", "gene_biotype"))
genes_names = genes_names[values(genes_names)[,"type"] == "gene"]
genes_names = genes_names[which(seqnames(genes_names) %in% c(1:20, "X")),]

genes_names_localQTL = genes_names[which(genes_names$gene_id %in% all_locally_regulated_genes),]

all_snps = read.table("/home/fkreuch/work/RI_translational_efficiency/data/20160318/RN6_SDPs/HXBBXH-SDP_RN60_markers.txt", sep="\t", comment.char = "#", header=F, stringsAsFactors=F, colClasses="character", check.names=F)

all_snps.GR = GRanges(seqnames = all_snps$V2, ranges = IRanges(start = as.numeric(all_snps$V3), end = as.numeric(all_snps$V3)), var = paste(all_snps$V4, all_snps$V5, sep="/"))

snps_in_locally_regulated_genes = subsetByOverlaps(all_snps.GR, genes_names_localQTL)

save.image("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/variant_effect_prediction/all_markers_locally_regulated_genes.RData")

write.table(data.frame(seqnames(snps_in_locally_regulated_genes), start(snps_in_locally_regulated_genes), start(snps_in_locally_regulated_genes), snps_in_locally_regulated_genes$var, 1), file = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/variant_effect_prediction/local_regGenes_input_VEP.txt", sep="\t", col.names=F, row.names=F, quote=F)




vep_genes = read.table("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/variant_effect_prediction/local_polyA_lv_VEP_output.txt", sep="\t", header=F, stringsAsFactors=F)



high_impact_lv_polyA = vep_lv_polyA[grep("HIGH", vep_lv_polyA$V5),]
high_impact_lv_ribo = vep_lv_ribo[grep("HIGH", vep_lv_ribo$V5),]
high_impact_lv_te = vep_lv_te[grep("HIGH", vep_lv_te$V5),]
high_impact_li_polyA = vep_li_polyA[grep("HIGH", vep_li_polyA$V5),]
high_impact_li_ribo = vep_li_ribo[grep("HIGH", vep_li_ribo$V5),]


load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/variant_effect_prediction/all_single_markers.RData")


lv_polyA = data.frame(sdp_id = lv_polyA_snp_ids$sdp, snp_id = lv_polyA_snp_ids$snp, conv = paste(all_snps[match(lv_polyA_snp_ids$snp, all_snps$V1), "V2"], all_snps[match(lv_polyA_snp_ids$snp, all_snps$V1), "V3"], paste(all_snps[match(lv_polyA_snp_ids$snp, all_snps$V1), "V4"], all_snps[match(lv_polyA_snp_ids$snp, all_snps$V1), "V5"],sep="/"), sep="_"))
lv_polyA = data.frame(lv_polyA, vep_lv_polyA[match(lv_polyA$conv, vep_lv_polyA$V1), c("V4", "V5", "V7", "V28")])
lv_polyA = data.frame(lv_polyA, lv_polyA_local[match(lv_polyA$sdp_id, lv_polyA_local$sdp_test),])
lv_polyA = lv_polyA[which(lv_polyA$V7 == lv_polyA$gene),]
#no high impact modifications in same gene as variant
#20 moderate impact -> 20 variants in 19 genes
lv_polyA = lv_polyA[which(lv_polyA$V5 == "MODERATE"),]
lv_polyA$QTL = "lv_polyA"


lv_ribo = data.frame(sdp_id = lv_ribo_snp_ids$sdp, snp_id = lv_ribo_snp_ids$snp, conv = paste(all_snps[match(lv_ribo_snp_ids$snp, all_snps$V1), "V2"], all_snps[match(lv_ribo_snp_ids$snp, all_snps$V1), "V3"], paste(all_snps[match(lv_ribo_snp_ids$snp, all_snps$V1), "V4"], all_snps[match(lv_ribo_snp_ids$snp, all_snps$V1), "V5"],sep="/"), sep="_"))
lv_ribo = data.frame(lv_ribo, vep_lv_ribo[match(lv_ribo$conv, vep_lv_ribo$V1), c("V4", "V5", "V7", "V28")])
lv_ribo = data.frame(lv_ribo, lv_ribo_local[match(lv_ribo$sdp_id, lv_ribo_local$sdp_test),])
lv_ribo = lv_ribo[which(lv_ribo$V7 == lv_ribo$gene),]
#no high impact modifications in same gene as variant
#14 moderate impact -> 14 variants in 11 genes
lv_ribo = lv_ribo[which(lv_ribo$V5 == "MODERATE"),]
lv_ribo$QTL = "lv_ribo"


lv_te = data.frame(sdp_id = lv_te_snp_ids$sdp, snp_id = lv_te_snp_ids$snp, conv = paste(all_snps[match(lv_te_snp_ids$snp, all_snps$V1), "V2"], all_snps[match(lv_te_snp_ids$snp, all_snps$V1), "V3"], paste(all_snps[match(lv_te_snp_ids$snp, all_snps$V1), "V4"], all_snps[match(lv_te_snp_ids$snp, all_snps$V1), "V5"],sep="/"), sep="_"))
lv_te = data.frame(lv_te, vep_lv_te[match(lv_te$conv, vep_lv_te$V1), c("V4", "V5", "V7", "V28")])
lv_te = data.frame(lv_te, lv_te_local[match(lv_te$sdp_id, lv_te_local$sdp_test),])
lv_te = lv_te[which(lv_te$V7 == lv_te$gene),]
#no high impact modifications in same gene as variant
#7 moderate impact -> 7 variants in 5 genes
lv_te = lv_te[which(lv_te$V5 == "MODERATE"),]
lv_te$QTL = "lv_te"




li_polyA = data.frame(sdp_id = li_polyA_snp_ids$sdp, snp_id = li_polyA_snp_ids$snp, conv = paste(all_snps[match(li_polyA_snp_ids$snp, all_snps$V1), "V2"], all_snps[match(li_polyA_snp_ids$snp, all_snps$V1), "V3"], paste(all_snps[match(li_polyA_snp_ids$snp, all_snps$V1), "V4"], all_snps[match(li_polyA_snp_ids$snp, all_snps$V1), "V5"],sep="/"), sep="_"))
li_polyA = data.frame(li_polyA, vep_li_polyA[match(li_polyA$conv, vep_li_polyA$V1), c("V4", "V5", "V7", "V28")])
li_polyA = data.frame(li_polyA, li_polyA_local[match(li_polyA$sdp_id, li_polyA_local$sdp_test),])
li_polyA = li_polyA[which(li_polyA$V7 == li_polyA$gene),]
#2 high impact modifications in same gene as variant - 1 splice_donor (Mst1), 1 stop gain (CES2e)
#11 moderate impact -> 11 variants in 11 genes
li_polyA = li_polyA[which(li_polyA$V5 %in% c("HIGH","MODERATE")),]
li_polyA$QTL = "li_polyA"


li_ribo = data.frame(sdp_id = li_ribo_snp_ids$sdp, snp_id = li_ribo_snp_ids$snp, conv = paste(all_snps[match(li_ribo_snp_ids$snp, all_snps$V1), "V2"], all_snps[match(li_ribo_snp_ids$snp, all_snps$V1), "V3"], paste(all_snps[match(li_ribo_snp_ids$snp, all_snps$V1), "V4"], all_snps[match(li_ribo_snp_ids$snp, all_snps$V1), "V5"],sep="/"), sep="_"))
li_ribo = data.frame(li_ribo, vep_li_ribo[match(li_ribo$conv, vep_li_ribo$V1), c("V4", "V5", "V7", "V28")])
li_ribo = data.frame(li_ribo, li_ribo_local[match(li_ribo$sdp_id, li_ribo_local$sdp_test),])
li_ribo = li_ribo[which(li_ribo$V7 == li_ribo$gene),]
#2 high impact modifications in same gene as variant, 
#8 moderate impact -> 8 variants in 8 genes
li_ribo = li_ribo[which(li_ribo$V5 %in% c("HIGH","MODERATE")),]
li_ribo$QTL = "li_ribo"

combine = rbind(lv_polyA,lv_ribo,lv_te,li_polyA,li_ribo)

write.table(combine, file = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/variant_effect_prediction/vep_localQTLs.txt", sep="\t", quote=F, col.names=T, row.names=F)

all_interesting_variants = data.frame(sdp = c("SDPG_19000024718", "SDPG_08114242338", "SDPG_20040480552", "SDPG_18072520144", "SDPG_02261053759", "SDPG_20012854759", "SDPG_10109484101", "SDPG_19041680327", "SDPG_02157592111"), gene = c("ENSRNOG00000011635", "ENSRNOG00000019680", "ENSRNOG00000000605", "ENSRNOG00000043171", "ENSRNOG00000028225", "ENSRNOG00000054549", "ENSRNOG00000036673", "ENSRNOG00000015063", "ENSRNOG00000012427"))

save(all_interesting_variants, file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/variant_effect_prediction/localQTLs_withImpact.RData")



