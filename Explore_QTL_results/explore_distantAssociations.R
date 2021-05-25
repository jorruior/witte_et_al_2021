
source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/SeqQTL_oldCluster.R")

library(GenomicRanges)
load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/permutations/All_tables_permutations.RData")


lv_polyA_distant

--> 3 genes regulated by SDPG_04017867424 (Dph5, Dram2, Myo5b)
--> 2 genes regulated by SDPG_07116912684 (Dusp23, Vamp7)
			SDPG_05130556071 (Sctr, Gli2)
			SDPG_02207507318 (Myo5b, Tmem158)
			SDPG_01254905321 (Tmcc2, Kel) ##Tmcc2 only gene that is distantly regulated in both tissues
			SDPG_01042979263 (Igfbp6, Atp5g2)

li_polyA_distant
--> 4 genes regulated by SDPG_06006563392 (Gas7, Ccdc64, Cyp2c13, Car12)
--> 2 genes regulated by SDPG_06006464104 (Ccdc64, S100a10)

lv_ribo_distant
--> 9 genes regulated by SDPG_03006314843 (Lama4, Lrp1,Notch2, Strn3, Lamb1, Fat4, Ssc5d, Agrn, AABR07073181.1)
--> 5 genes regulated by SDPG_03003960268 (Lama4, Lrp1, Fat4, Notch2, Lamb1)
--> 4 genes regulated by SDPG_08051402844 (Rasgrp1, Atp6ap2, Strbp, Dlgap1)
--> 3 genes regulated by SDPG_08052253560 (Nr2c2ap, Ost4, Phf3)
			SDPG_03009805651 (Fat4, Lrp1, Lama4)
--> 2 genes regulated by SDPG_08095189269 (Emp1, Bri3bp)
			SDPG_08065229351 (Mlec, Dtnbp1)			

li_ribo_distant
--> no associations

lv_te_distant
--> 11 genes regulated by SDPG_03003960268 (Agrn, Lrp1, Dchs1, Kdr, Hmcn1, Lam4, Igf2r, Fat1, Lamb1, Sema3c, Uggt1)
--> 10 genes regulated by SDPG_03006314843 (Agrn, Lrp1, Dchs1,AABR07073181.1, Lam4, Fat1, Igf2r, Lamb1, Sema3c, Stt3a) 
--> 8 genes regulated by SDPG_08051402844 (Enpp5, Atp2b2, Acot7, Tspan7, Fam171b, Trim37, Ppp1r9a, Nlgn2)
--> 6 genes regulated by SDPG_03009805651 (Kdr, Lrp1, Hmcn1, Dchs1, Agrn, Igf2r)
--> 4 genes regulated by SDPG_08059217527 (Rnpepl1, Mfhas1, Nolc1, Eif4g3)
			SDPG_03000054231 (Stt3a, Sema3c, LOC679811, Hmcn1)
--> 3 genes regulated by SDPG_10106181521 (Itgb1, Heph, B4galt1)
			SDPG_08052253560 (Znrf1, Ost4, Cmss1)
			SDPG_05079069190 (Ano5, Arhgab26, Fam220a)
--> 2 genes regulated by SDPG_11039674052 (Mmd, Ddx27)
			SDPG_08066825758 (Tex261, Wdr45)
			SDPG_08066150922 (Scarf1, Klhl13)
			SDPG_08065114132 (Ube2q2l, Wdr45)
			SDPG_04161793548 (Trem2, Fjx1)
			SDPG_02045317417 (Tmem150a, Shisa5)
			SDPG_02044100706 (Tmem150a, Shisa5)






genotypes = read.table("~/work/RI_translational_efficiency/data/20160318/RN6_SDPs/HXBBXH-SDP_RN60_LDB-genome.txt",  sep="\t", comment.char = "#", header=F, stringsAsFactors=F, colClasses="character", check.names=F)
colnames(genotypes) = c("ID", "RI_geno", "SDP_local")

####kick out single-marker SDPs
local_sdps = read.table("/home/fkreuch/work/RI_translational_efficiency/data/20160318/RN6_SDPs/HXBBXH-SDP_RN60_LDB-local.txt", sep="\t", comment.char = "#", header=F, stringsAsFactors=F, colClasses="character", check.names=F)
colnames(local_sdps) = c("local_ID", "chr", "start", "end", "Marker_RN60_CSV", "Marker_count")


local_sdps[c(which(local_sdps$local_ID == "SDPL_08051402844"):which(local_sdps$local_ID == "SDPL_08095189269")), c(1:4,6) ]

sample_names = c("BXH02", "BXH03", "BXH05", "HXB01", "HXB02", "HXB03", "HXB04", "HXB05", "HXB07", "HXB10", "HXB13", "HXB15", "HXB17", "HXB18", "HXB20", "HXB21", "HXB22", "HXB23", "HXB24", "HXB25", "HXB27", "HXB29", "HXB31", "BXH06", "BXH08", "BXH09", "BXH10", "BXH11", "BXH12", "BXH13") 
chr8_genotypePatterns = data.frame(apply(genotypes[c(which(genotypes$ID == "SDPG_08051402844") : which(genotypes$ID == "SDPG_08095189269")),], 1, function(x){strsplit(x[2], split="")}))
genotypes_alleles = t(chr8_genotypePatterns[c(7:36),])
colnames(genotypes_alleles) = sample_names
rownames(genotypes_alleles) = genotypes[c(which(genotypes$ID == "SDPG_08051402844") : which(genotypes$ID == "SDPG_08095189269")),"ID"]


library(ggplot2)
library(reshape2)
melt_chr8_geno = melt(genotypes_alleles)
pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/plots/chr8_clusterSDPs.pdf", width=10, height=7)
ggplot(data= melt_chr8_geno , aes(x=Var2, y=Var1, col = factor(value)))+geom_point(aes(size=2))+theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

















all_sets = list(lv_polyA_local = unique(lv_polyA_local$gene), lv_polyA_distant = unique(lv_polyA_distant$gene), lv_ribo_local = unique(lv_ribo_local$gene), lv_ribo_distant = unique(lv_ribo_distant$gene), lv_te_distant = unique(lv_te_distant$gene), lv_te_local = unique(lv_te_local$gene), li_polyA_local = unique(li_polyA_local$gene), li_polyA_distant = unique(li_polyA_distant$gene), li_ribo_local = unique(li_ribo_local$gene), li_ribo_distant = unique(li_ribo_distant$gene), li_te_distant = unique(li_te_distant$gene), li_te_local = unique(li_te_local$gene))


#same results
library(gplots)
ItemsList <- venn(all_sets, show.plot=FALSE)
intersections = intersections<-attr(ItemsList,"intersections")



ItemsList_dist_lv <- venn(list(rna = unique(lv_polyA_distant$gene), ribo = unique(lv_ribo_distant$gene), te = unique(lv_te_distant$gene) ), show.plot=T)



lv_polyA_cluster = lv_polyA_distant[which(lv_polyA_distant$SNP %in% names(table(lv_polyA_distant$SNP))[table(lv_polyA_distant$SNP) >=2]), c("SNP", "gene", "gene_name")]
li_polyA_cluster = li_polyA_distant[which(li_polyA_distant$SNP %in% names(table(li_polyA_distant$SNP))[table(li_polyA_distant$SNP) >=2]), c("SNP", "gene", "gene_name")]

lv_ribo_cluster = lv_ribo_distant[which(lv_ribo_distant$SNP %in% names(table(lv_ribo_distant$SNP))[table(lv_ribo_distant$SNP) >=2]), c("SNP", "gene", "gene_name")]
li_ribo_cluster = li_ribo_distant[which(li_ribo_distant$SNP %in% names(table(li_ribo_distant$SNP))[table(li_ribo_distant$SNP) >=2]), c("SNP", "gene", "gene_name")]

lv_te_cluster = lv_te_distant[which(lv_te_distant$SNP %in% names(table(lv_te_distant$SNP))[table(lv_te_distant$SNP) >=2]), c("SNP", "gene", "gene_name")]
li_te_cluster = li_te_distant[which(li_te_distant$SNP %in% names(table(li_te_distant$SNP))[table(li_te_distant$SNP) >=2]), c("SNP", "gene", "gene_name")]


all_poss_cluster = rbind(lv_polyA_cluster, li_polyA_cluster, lv_ribo_cluster, li_ribo_cluster, lv_te_cluster, li_te_cluster)
all_poss_cluster$asso = as.character(sapply(all_poss_cluster$gene, function(x) names(intersections)[grep(x ,intersections)]))
write.table(all_poss_cluster, file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/distantCluster/distant_cluster.txt", sep="\t", col.names=T, row.names=F, quote=F)






load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/permutations/All_tables_permutations.RData")
load(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/preprocessed_data.RData")


sample_names = c("BXH02", "BXH03", "BXH05", "HXB01", "HXB02", "HXB03", "HXB04", "HXB05", "HXB07", "HXB10", "HXB13", "HXB15", "HXB17", "HXB18", "HXB20", "HXB21", "HXB22", "HXB23", "HXB24", "HXB25", "HXB27", "HXB29", "HXB31", "BXH06", "BXH08", "BXH09", "BXH10", "BXH11", "BXH12", "BXH13") 

#Phenotype Input file for eQTL Mapping
log2_polyA_lv = log2_polyA_lv[,match(sample_names, colnames(log2_polyA_lv))]
log2_polyA_li = log2_polyA_li[,match(sample_names, colnames(log2_polyA_li))]
log2_ribo_lv = log2_ribo_lv[,match(sample_names, colnames(log2_ribo_lv))]
colnames(log2_ribo_li) = gsub("_4", "", colnames(log2_ribo_li))
log2_ribo_li = log2_ribo_li[,match(sample_names, colnames(log2_ribo_li))]

log2_te_lv = resid(lm(as.matrix(log2_ribo_lv)~as.matrix(log2_polyA_lv)))
log2_te_li = resid(lm(as.matrix(log2_ribo_li)~as.matrix(log2_polyA_li)))

distantCluster = read.table(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/distantCluster/distant_cluster.txt", sep="\t", header=T, stringsAsFactors=F)

snp = c("SDPG_03003960268", "SDPG_03006314843", "SDPG_03009805651")
genes_of_Interest = unique(distantCluster[which(distantCluster$SNP %in% snp), "gene"])

correlations_cluster2 = cor(t(log2_te_lv), t(log2_te_lv[genes_of_Interest,]), method="spearman")
genes_names = gtf2GR("/data/huebner2/DATASETS/GENOMES/Ensembl_Genomes/Ens82/rn6/annotation/Rattus_norvegicus.Rnor_6.0.82.gtf", c("gene_id", "gene_name", "gene_biotype"))
genes_names = genes_names[values(genes_names)[,"type"] == "gene"]
genes_names = genes_names[which(seqnames(genes_names) %in% c(1:20, "X")),]

rownames(correlations_cluster2) = genes_names[match(rownames(correlations_cluster2), genes_names$gene_id), ]$gene_name
colnames(correlations_cluster2) = genes_names[match(colnames(correlations_cluster2), genes_names$gene_id), ]$gene_name


library(rtracklayer)
library('dendextend')
library(gplots)
library(HTSCluster)
library(fields)

hc_qtl <- hclust(dist(t(abs(correlations_cluster2)))) 
hc_qtl <-as.dendrogram(hc_qtl)
hc_qtl %>%  plot(main = "Change colors") 
ct_qtl <- cutree(hc_qtl, k=3) 
qtl_cut=rect.hclust(as.hclust(hc_qtl), k=3,border = ct_qtl) 


hc_all <- hclust(dist(abs(correlations_cluster2))) 
hc_all <-as.dendrogram(hc_all)
hc_all %>%  plot(main = "Change colors") 
ct_all <- cutree(hc_all, k=7) 
all_cut=rect.hclust(as.hclust(hc_all), k=7,border = ct_all) 


colors = rainbow(16, s = 0.5)
mycols  <- colorRampPalette(c("blue", "white", "red"), space = "Lab")(n = 20)

pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/distantCluster/test_correlation_17Cluster2Genevsall.pdf", width=12,height=12 )


heatmap.2(correlations_cluster2,trace='none',Rowv=as.dendrogram(hc_all), Colv=as.dendrogram(hc_qtl), density.info='none',
             col=mycols,key = TRUE,ColSideColors = colors[as.factor(ct_qtl)],RowSideColors = colors[as.factor(ct_all)],
             scale="none",labCol=colnames(correlations_cluster2),labRow=NA,srtCol=45)



dev.off()


#cluster 6 and 2
subset_genes_high_correlation = names(ct_all[ct_all %in% c(2,6)])
subset_correlations_cluster2 = correlations_cluster2[match(subset_genes_high_correlation, rownames(correlations_cluster2)),]
hc_qtl2 <- hclust(dist(t(abs(subset_correlations_cluster2)))) 
hc_qtl2 <-as.dendrogram(hc_qtl2)
hc_qtl2 %>%  plot(main = "Change colors") 
ct_qtl2 <- cutree(hc_qtl2, k=3) 
qtl_cut2=rect.hclust(as.hclust(hc_qtl2), k=3,border = ct_qtl2) 


hc_all2 <- hclust(dist(abs(subset_correlations_cluster2))) 
hc_all2 <-as.dendrogram(hc_all2)
hc_all2 %>%  plot(main = "Change colors") 
ct_all2 <- cutree(hc_all2, k=7) 
all_cut2=rect.hclust(as.hclust(hc_all2), k=7,border = ct_all2) 


colors = rainbow(16, s = 0.5)
mycols  <- colorRampPalette(c("blue", "white", "red"), space = "Lab")(n = 20)

pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/distantCluster/test_correlation_17Cluster2GenevsSubset.pdf", width=12,height=12 )


heatmap.2(subset_correlations_cluster2,trace='none',Rowv=as.dendrogram(hc_all2), Colv=as.dendrogram(hc_qtl2), density.info='none',
             col=mycols,key = TRUE,ColSideColors = colors[as.factor(ct_qtl2)],RowSideColors = colors[as.factor(ct_all2)],
             scale="none",labCol=colnames(subset_correlations_cluster2),labRow=rownames(subset_correlations_cluster2),srtCol=45)

dev.off()


scaled_counts = t(scale(t(log2_te_lv)))

genotypes_alleles = read.table(file = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Input/sdp_genotypes.txt", sep="\t", header=T, stringsAsFactors=F)
snp = "SDPG_03003960268"
pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/distantCluster/plots_clustered_groups_of_expression.pdf")

cluster_genes = genes_names[match(rownames(subset_correlations_cluster2), genes_names$gene_name),]$gene_id

genes_toPlot_cluster = as.data.frame(scaled_counts[which(rownames(scaled_counts) %in% cluster_genes),])
genes_toPlot_cluster$gene_id = rownames(genes_toPlot_cluster )
genes_toPlot_cluster$group = "gene"
geno = genotypes_alleles[snp,]

genes_toPlot_cluster[which(genes_toPlot_cluster$gene_id %in% genes_of_Interest),]$group = "QTL_genes"


library(reshape)
genes_toPlot_cluster.melted = melt(genes_toPlot_cluster, id = c("gene_id", "group"))
genes_toPlot_cluster.melted$genotype = as.character(t(geno[,match(genes_toPlot_cluster.melted$variable, colnames(geno))]))
library(ggplot2)


test = ggplot(data=genes_toPlot_cluster.melted[which(genes_toPlot_cluster.melted$group == "gene"),], aes(x=variable, y=value, group = gene_id))+geom_point(color="lightgrey")+ylab("")+ggtitle("Scaled Expression Cluster")+facet_grid(.~genotype)

test = test + geom_point(aes(variable, value, group=gene_id), color = "darksalmon", genes_toPlot_cluster.melted[which(genes_toPlot_cluster.melted$group == "QTL_genes"),])

print(test)

}
dev.off()







subset_correlations = sapply(colnames(correlations_cluster2), function(x) correlations_cluster2[which(abs(correlations_cluster2[,x]) > 0.8), x])

load(file = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/mirBase/rat_mirBindingTargets.RData")
cluster2_genes_with_miRNA_site = mirbase_rat.min[which(mirbase_rat.min$gene_id %in% genes_of_Interest),]
cluster2_genes_with_miRNA_site$gene_name = genes_names[match(cluster2_genes_with_miRNA_site$gene_id, genes_names$gene_id), ]$gene_name
#nothing enriched/common


mtor="ENSRNOG00000009615"




