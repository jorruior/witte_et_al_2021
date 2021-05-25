load.counts_CDS <- function(files, ann) {
  counts = NULL
  for (file in files) {
    x = read.table(file, header=F, stringsAsFactors=F)
    x = data.frame(id = x$V1, x$V2)

    colnames(x) = c("id", gsub("_CDS_counts.txt", "", basename(file)))
    if (is.null(counts)) {
      counts = x
    } else {
      counts = merge(counts, x, all.x=T, all.y=T)
    }
  }
  rownames(counts) = counts[,"id"]
  counts = counts[,-1]
   
  missing.exons = setdiff(elementMetadata(ann)[,"gene_id"], rownames(counts))
  missing.counts = matrix(ncol=ncol(counts), nrow=length(missing.exons))
  colnames(missing.counts) = colnames(counts)
  rownames(missing.counts) = missing.exons
  
  counts = rbind(counts, missing.counts)
  counts[is.na(counts)] = 0
  counts = counts[elementMetadata(ann)[,"gene_id"],]
  return(counts)
}
library(rtracklayer)
cols=c("type","source","gene_id","gene_version","transcript_id","exon_number","gene_name","gene_source","gene_biotype","havana_gene",
			"havana_gene_version","transcript_name","transcript_source","transcript_biotype","exon_id","exon_version")
library(ggplot2)
library(beeswarm)


mORF_uORF_ranges_lv = import.gff("/Volumes//huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/uORFs_quant_complete_lv.gtf", colnames=cols)
overlapping_uORF_ranges_lv = import.gff("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/uORFs_quant_only_overlapping_uORFs_lv.gtf", colnames=cols)

mORF_uORF_ranges_li = import.gff("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/uORFs_quant_complete_li.gtf", colnames=cols)
overlapping_uORF_ranges_li = import.gff("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/uORFs_quant_only_overlapping_uORFs_li.gtf", colnames=cols)


mORF_uORF_ranges_names_lv = subset(mORF_uORF_ranges_lv, !duplicated(mORF_uORF_ranges_lv$gene_id))
overlapping_uORF_ranges_names_lv = subset(overlapping_uORF_ranges_lv, !duplicated(overlapping_uORF_ranges_lv$gene_id))
mORF_uORF_ranges_names_li = subset(mORF_uORF_ranges_li, !duplicated(mORF_uORF_ranges_li$gene_id))
overlapping_uORF_ranges_names_li = subset(overlapping_uORF_ranges_li, !duplicated(overlapping_uORF_ranges_li$gene_id))

files_ribo_CDS_lv_mORF_uORF = list.files("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/lv_ribo_indep_uORF" , pattern= "_CDS_counts.txt", full.names=TRUE)
files_ribo_CDS_lv_overlapping_uORF = list.files("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/lv_ribo" , pattern= "_CDS_counts.txt", full.names=TRUE)
files_ribo_CDS_li_mORF_uORF = list.files("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/li_ribo_indep_uORF" , pattern= "_CDS_counts.txt", full.names=TRUE)
files_ribo_CDS_li_overlapping_uORF = list.files("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/li_ribo" , pattern= "_CDS_counts.txt", full.names=TRUE)


counts_ribo_CDS_overlapping_uORF_lv = load.counts_CDS(files_ribo_CDS_lv_overlapping_uORF, overlapping_uORF_ranges_names_lv)
counts_ribo_CDS_mORF_uORF_lv = load.counts_CDS(files_ribo_CDS_lv_mORF_uORF, mORF_uORF_ranges_names_lv)
counts_ribo_CDS_overlapping_uORF_li = load.counts_CDS(files_ribo_CDS_li_overlapping_uORF, overlapping_uORF_ranges_names_li)
counts_ribo_CDS_mORF_uORF_li = load.counts_CDS(files_ribo_CDS_li_mORF_uORF, mORF_uORF_ranges_names_li)

load("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_tables.RData")

filtered_genes_lv = unique(all_uORFs_lv$gene_id)
filtered_genes_li = unique(all_uORFs_li$gene_id)

load(file="/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/filtered_genesets_ribotaper_fpkm.RData")

all_uORFs_li = all_uORFs_li[which(all_uORFs_li$category %in% c("uORF", "ncORFS", "ORFs_ccds")),]

filtered_genes_lv = filtered_genes_lv[which(filtered_genes_lv %in% genesset_lv)]
filtered_genes_li = filtered_genes_li[which(filtered_genes_li %in% genesset_li)]

all_uORFs_lv.filt = all_uORFs_lv[which(all_uORFs_lv$gene_id %in% filtered_genes_lv),]
all_uORFs_li.filt = all_uORFs_li[which(all_uORFs_li$gene_id %in% filtered_genes_li),]

length_normalization = function(counts, annotation, total)
{
bygene = split(annotation, values(annotation)[, "gene_id"])
len = sum(width(reduce(bygene)))
genes = intersect(names(bygene), rownames(counts))
counts = counts[genes, ]
len = len[genes]

len_norm_counts = counts/rep(total, each = nrow(counts)) / rep(len, ncol(counts)) * 1e+09

return(len_norm_counts)

}

ribo_lv_counts_raw = read.table(file="/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/raw/ribo_RI_lv_proteincodingANDlincRNA_countMatrix_raw.txt", sep="\t", header=T, stringsAsFactors=F)
number_total_counts_lv = colSums(ribo_lv_counts_raw)

ribo_li_counts_raw = read.table(file="/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/raw/ribo_RI_li_proteincodingANDlincRNA_countMatrix_raw.txt", sep="\t", header=T, stringsAsFactors=F)
number_total_counts_li = colSums(ribo_li_counts_raw)


counts_ribo_CDS_uORF_nonoverlapping_lv = counts_ribo_CDS_mORF_uORF_lv[grep("independent", rownames(counts_ribo_CDS_mORF_uORF_lv)),]

lenNorm_ribo_uORF_nonOverlapping_lv = length_normalization(counts_ribo_CDS_uORF_nonoverlapping_lv, mORF_uORF_ranges_lv, number_total_counts_lv)
lenNorm_ribo_uORF_overlap_lv = length_normalization(counts_ribo_CDS_overlapping_uORF_lv, overlapping_uORF_ranges_lv, number_total_counts_lv)

counts_ribo_CDS_uORF_nonoverlapping_li = counts_ribo_CDS_mORF_uORF_li[grep("independent", rownames(counts_ribo_CDS_mORF_uORF_li)),]
lenNorm_ribo_uORF_nonOverlapping_li = length_normalization(counts_ribo_CDS_uORF_nonoverlapping_li, mORF_uORF_ranges_li, number_total_counts_li)
lenNorm_ribo_uORF_overlap_li = length_normalization(counts_ribo_CDS_overlapping_uORF_li, overlapping_uORF_ranges_li, number_total_counts_li)

lenNorm_ribo_uORF_overlap_li$ORF_id_tr = overlapping_uORF_ranges_names_li[match(rownames(lenNorm_ribo_uORF_overlap_li),overlapping_uORF_ranges_names_li$gene_id ), ]$transcript_id
lenNorm_ribo_uORF_overlap_lv$ORF_id_tr = overlapping_uORF_ranges_names_lv[match(rownames(lenNorm_ribo_uORF_overlap_lv),overlapping_uORF_ranges_names_lv$gene_id ), ]$transcript_id

lenNorm_ribo_uORF_nonOverlapping_li$ORF_id_tr = mORF_uORF_ranges_names_li[match(rownames(lenNorm_ribo_uORF_nonOverlapping_li),mORF_uORF_ranges_names_li$gene_id ), ]$transcript_id
lenNorm_ribo_uORF_nonOverlapping_lv$ORF_id_tr = mORF_uORF_ranges_names_lv[match(rownames(lenNorm_ribo_uORF_nonOverlapping_lv),mORF_uORF_ranges_names_lv$gene_id ), ]$transcript_id


length_Norm_ribo_uORF_lv = rbind(lenNorm_ribo_uORF_nonOverlapping_lv,lenNorm_ribo_uORF_overlap_lv)
length_Norm_ribo_uORF_li = rbind(lenNorm_ribo_uORF_nonOverlapping_li,lenNorm_ribo_uORF_overlap_li)


sample_names = c("BXH02", "BXH03", "BXH05", "HXB01", "HXB02", "HXB03", "HXB04", "HXB05", "HXB07", "HXB10", "HXB13", "HXB15", "HXB17", "HXB18", "HXB20", "HXB21", "HXB22", "HXB23", "HXB24", "HXB25", "HXB27", "HXB29", "HXB31", "BXH06", "BXH08", "BXH09", "BXH10", "BXH11", "BXH12", "BXH13") 

length_Norm_ribo_uORF_lv = length_Norm_ribo_uORF_lv[,match(c(sample_names, "ORF_id_tr"), gsub("_lv", "", colnames(length_Norm_ribo_uORF_lv)))]

length_Norm_ribo_uORF_li = length_Norm_ribo_uORF_li[,match(c(sample_names, "ORF_id_tr"), gsub("_4", "",gsub("_li", "", colnames(length_Norm_ribo_uORF_li))))]


load(file = "/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Input/Lv_Ribo_RNA_resid/Lv_Ribo_RNA_resid_input.RData")
TE_lv = pheno
load(file = "/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Input/Liver_Ribo_RNA_resid/Liver_Ribo_RNA_resid_input.RData")
TE_li = pheno


library(reshape2)

var_TE_lv = apply(TE_lv, 1, var)
melted_TE_lv = melt(TE_lv)
melted_TE_lv$category = "no uORF"
melted_TE_lv[which(melted_TE_lv$Var1 %in% all_uORFs_lv[which(all_uORFs_lv$category == "uORF"),"gene_id"]),]$category = "independent uORF"
melted_TE_lv[which(melted_TE_lv$Var1 %in% all_uORFs_lv[which(all_uORFs_lv$category != "uORF"),"gene_id"]),]$category = "overlapping uORF"

variance_lv = data.frame(gene_id = names(var_TE_lv), value = var_TE_lv)
variance_lv$category = "no uORF"
variance_lv[which(variance_lv$gene_id %in% all_uORFs_lv[which(all_uORFs_lv$category == "uORF"),"gene_id"]),]$category = "independent uORF"
variance_lv[which(variance_lv$gene_id %in% all_uORFs_lv[which(all_uORFs_lv$category != "uORF"),"gene_id"]),]$category = "overlapping uORF"

var_TE_li = apply(TE_li, 1, var)
melted_TE_li = melt(TE_li)
melted_TE_li$category = "no uORF"
melted_TE_li[which(melted_TE_li$Var1 %in% all_uORFs_li[which(all_uORFs_li$category == "uORF"),"gene_id"]),]$category = "independent uORF"
melted_TE_li[which(melted_TE_li$Var1 %in% all_uORFs_li[which(all_uORFs_li$category != "uORF"),"gene_id"]),]$category = "overlapping uORF"

variance_li = data.frame(gene_id = names(var_TE_li), value = var_TE_li)
variance_li$category = "no uORF"
variance_li[which(variance_li$gene_id %in% all_uORFs_li[which(all_uORFs_li$category == "uORF"),"gene_id"]),]$category = "independent uORF"
variance_li[which(variance_li$gene_id %in% all_uORFs_li[which(all_uORFs_li$category != "uORF"),"gene_id"]),]$category = "overlapping uORF"

library(ggplot2)
pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment/TE_mORF_impact_uORF_presence.pdf")
ggplot(melted_TE_lv, aes(x=as.factor(category), y=log2(value), fill=as.factor(category)))+geom_boxplot(width=0.5)+ ggtitle("Impact of uORF on TE of mORF LV")+ylab("log2(TE)")+xlab("")
ggplot(variance_lv, aes(x=as.factor(category), y=log2(value), fill=as.factor(category)))+geom_boxplot(width=0.5)+ ggtitle("Impact of uORF on variance of TE of mORF LV")+ylab("log2(variance TE)")+xlab("")
ggplot(melted_TE_li, aes(x=as.factor(category), y=log2(value), fill=as.factor(category)))+geom_boxplot(width=0.5)+ ggtitle("Impact of uORF on TE of mORF Liver")+ylab("log2(TE)")+xlab("")
ggplot(variance_li, aes(x=as.factor(category), y=log2(value), fill=as.factor(category)))+geom_boxplot(width=0.5)+ ggtitle("Impact of uORF on variance of TE of mORF Liver")+ylab("log2(variance TE)")+xlab("")
dev.off()


save.image("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment/20180227_image.RData")

length_norm_uORF_lv = data.frame(length_Norm_ribo_uORF_lv, gene_id = gsub("_overlapping", "", gsub("_independent", "", rownames(length_Norm_ribo_uORF_lv))))
length_norm_uORF_li = data.frame(length_Norm_ribo_uORF_li, gene_id = gsub("_overlapping", "", gsub("_independent", "", rownames(length_Norm_ribo_uORF_li))))

TE_of_uORF_corresp_genes_lv = TE_lv[match(as.character(length_norm_uORF_lv$gene_id), rownames(TE_lv)),]
TE_of_uORF_corresp_genes_li = TE_li[match(as.character(length_norm_uORF_li$gene_id), rownames(TE_li)),]

cor_mORF_TE_lengthNorm_uORF_ribo_lv = sapply(1:nrow(length_norm_uORF_lv), function(x) cor(t(length_norm_uORF_lv[x,1:30]), TE_of_uORF_corresp_genes_lv[x,], method="spearman"))

cor_mORF_TE_lengthNorm_uORF_ribo_li = sapply(1:nrow(length_norm_uORF_li), function(x) cor(t(length_norm_uORF_li[x,1:30]), TE_of_uORF_corresp_genes_li[x,], method="spearman"))

correlation_TE_mORF_length_uORF_lv = data.frame(gene_id = rownames(TE_of_uORF_corresp_genes_lv), value = cor_mORF_TE_lengthNorm_uORF_ribo_lv)
correlation_TE_mORF_length_uORF_lv$category = "independent"
correlation_TE_mORF_length_uORF_lv[which(correlation_TE_mORF_length_uORF_lv$gene_id %in% all_uORFs_lv[which(all_uORFs_lv$category != "uORF"),"gene_id"]),]$category = "overlapping"


correlation_TE_mORF_length_uORF_li = data.frame(gene_id = rownames(TE_of_uORF_corresp_genes_li), value = cor_mORF_TE_lengthNorm_uORF_ribo_li)
correlation_TE_mORF_length_uORF_li$category = "independent"
correlation_TE_mORF_length_uORF_li[which(correlation_TE_mORF_length_uORF_li$gene_id %in% all_uORFs_li[which(all_uORFs_li$category != "uORF"),"gene_id"]),]$category = "overlapping"


correlation_TE_mORF_length_uORF_lv$gene_name = all_uORFs_lv[match(correlation_TE_mORF_length_uORF_lv$gene_id, all_uORFs_lv$gene_id), "gene_symbol"]
correlation_TE_mORF_length_uORF_li$gene_name = all_uORFs_li[match(correlation_TE_mORF_length_uORF_li$gene_id, all_uORFs_li$gene_id), "gene_symbol"]


all_uORFs_lv = data.frame(all_uORFs_lv, correlation_TE_mORF_length_uORF_lv[match(all_uORFs_lv$gene_id, correlation_TE_mORF_length_uORF_lv$gene_id),])
all_uORFs_li = data.frame(all_uORFs_li, correlation_TE_mORF_length_uORF_li[match(all_uORFs_li$gene_id, correlation_TE_mORF_length_uORF_li$gene_id),])
write.table(all_uORFs_lv, file = "Desktop/Manuscript/uORF_table_left_ventricle.txt", sep="\t", col.names=T, row.names=F, quote = F)
write.table(all_uORFs_li, file = "Desktop/Manuscript/uORF_table_liver.txt", sep="\t", col.names=T, row.names=F, quote = F)

#> correlation_TE_mORF_length_uORF_lv[which(correlation_TE_mORF_length_uORF_lv$value < -0.5),]
               gene_id      value    category gene_name
28  ENSRNOG00000001114 -0.5323693 independent     Wipi2
59  ENSRNOG00000002635 -0.5817575 independent      Dexi
104 ENSRNOG00000004827 -0.5305895 independent    Papola
142 ENSRNOG00000006833 -0.5657397 independent    Rb1cc1
155 ENSRNOG00000007681 -0.5261402 independent      Brd3
299 ENSRNOG00000012844 -0.5959956 independent      Tox4
534 ENSRNOG00000024886 -0.5083426 overlapping      Ext1
620 ENSRNOG00000042492 -0.5692992 independent     Rfwd2

#> correlation_TE_mORF_length_uORF_li[which(correlation_TE_mORF_length_uORF_li$value < -0.5),]
               gene_id      value    category gene_name
86  ENSRNOG00000003144 -0.5328142 independent    Gprc5c
151 ENSRNOG00000005689 -0.5430478 independent    Yeats4
610 ENSRNOG00000021726 -0.5190211 independent      Tlr3
648 ENSRNOG00000025143 -0.6441206 independent     Icam2









correlation_TE_mORF_length_uORF_lv$name =  "Spearman rho"
correlation_TE_mORF_length_uORF_li$name =  "Spearman rho"

pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment//correlation_TE_uORF_plotsdiffGenes.pdf", width=10)

par(mfrow=c(1,3))

beeswarm(value ~ name, data=correlation_TE_mORF_length_uORF_lv, ylim =c(-1,1), pch=20, main = "Correlation mORF TE and uORF Ribo LV")
abline(h = 0)
beeswarm(value ~ name, data=correlation_TE_mORF_length_uORF_li, ylim =c(-1,1), pch=20, main = "Correlation mORF TE and uORF Ribo Liver")
abline(h = 0)

dev.off()



