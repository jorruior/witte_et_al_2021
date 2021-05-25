library(stringr)

load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/permutations/All_tables_permutations.RData")
load(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/filtered_genesets_ribotaper_fpkm.RData")


expressed_in_both = intersect(genesset_li, genesset_lv) #8622 genes
lv_specific = setdiff(genesset_lv, genesset_li) # 1909
li_specific = setdiff(genesset_li, genesset_lv) # 714

library(gplots)



#local QTLs

eQTL_local_same_gene = intersect(lv_polyA_local$gene, li_polyA_local$gene) #128 genes
eQTL_local_same_SDP_gene = intersect(paste(lv_polyA_local$gene, lv_polyA_local$SNP, sep="_"), paste(li_polyA_local$gene, li_polyA_local$SNP, sep="_")) #291 pairs
length(unique(str_split_fixed(eQTL_local_same_SDP_gene, "_", 3)[,1])) 

library(gProfileR)
set_base_url("http://biit.cs.ut.ee/gprofiler_archive2/r1477_e82_eg29/web")

polyA_shared = gprofiler(query = unique(eQTL_local_same_gene), custom_bg = expressed_in_both, organism = "rnorvegicus")

TF = read.table("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/transciption_factor_binding.txt", sep="", header=F, stringsAsFactors=F)


# 126 genes regulated by same SDP in both tissue --> 1.46%

eQTL_local_specific_lv = setdiff(lv_polyA_local$gene, li_polyA_local$gene) # 574 genes
eQTL_local_specific_li = setdiff(li_polyA_local$gene, lv_polyA_local$gene) # 237 genes

eQTL_local_only_expressed_lv = intersect(lv_polyA_local$gene, lv_specific) #142 genes with QTL only expressed in lv --> 7.44%
eQTL_local_only_expressed_li = intersect(li_polyA_local$gene, li_specific) #69 genes with QTL only expressed in li --> 8.26%

riboQTL_local_same_gene = intersect(lv_ribo_local$gene, li_ribo_local$gene) #85 genes
riboQTL_local_same_SDP_gene = intersect(paste(lv_ribo_local$gene, lv_ribo_local$SNP, sep="_"), paste(li_ribo_local$gene, li_ribo_local$SNP, sep="_")) #5 pairs
length(unique(str_split_fixed(riboQTL_local_same_SDP_gene, "_", 3)[,1])) # 81 genes regulated by same SDP in both tissue --> 0.94
riboQTL_local_specific_lv = setdiff(lv_ribo_local$gene, li_ribo_local$gene) # 224 genes
riboQTL_local_specific_li = setdiff(li_ribo_local$gene, lv_ribo_local$gene) # 298 genes

riboQTL_local_only_expressed_lv = intersect(lv_ribo_local$gene, lv_specific) #82 genes with QTL only expressed in lv --> 4.30%
riboQTL_local_only_expressed_li = intersect(li_ribo_local$gene, li_specific) #62 genes with QTL only expressed in li --> 8.68%


ribo_shared = gprofiler(query = unique(riboQTL_local_same_gene), custom_bg = expressed_in_both, organism = "rnorvegicus")


rbp = read.table("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/rbp.txt", sep="", header=F, stringsAsFactors=F)


teQTL_local_same_gene = intersect(lv_te_local$gene, li_te_local$gene) #22 genes
teQTL_local_same_SDP_gene = intersect(paste(lv_te_local$gene, lv_te_local$SNP, sep="_"), paste(li_te_local$gene, li_te_local$SNP, sep="_")) #5 pairs
length(unique(str_split_fixed(teQTL_local_same_SDP_gene, "_", 3)[,1])) # 21 genes regulated by same SDP in both tissue --> 0.24%

teQTL_local_specific_lv = setdiff(lv_te_local$gene, li_te_local$gene) # 49 genes
teQTL_local_specific_li = setdiff(li_te_local$gene, lv_te_local$gene) # 49 genes

teQTL_local_only_expressed_lv = intersect(lv_te_local$gene, lv_specific) #11 genes with QTL only expressed in lv --> 0.58%
teQTL_local_only_expressed_li = intersect(li_te_local$gene, li_specific) #0 genes with QTL only expressed in li --> 0%



ItemsList_lv <- venn(list(polyA = unique(str_split_fixed(eQTL_local_same_SDP_gene, "_", 3)[,1]),
ribo = unique(str_split_fixed(riboQTL_local_same_SDP_gene, "_", 3)[,1]),
te = unique(str_split_fixed(teQTL_local_same_SDP_gene, "_", 3)[,1])), show.plot=T)

#distant QTLs

eQTL_distant_same_gene = intersect(lv_polyA_distant$gene, li_polyA_distant$gene) #1 genes
eQTL_distant_same_SDP_gene = intersect(paste(lv_polyA_distant$gene, lv_polyA_distant$SNP, sep="_"), paste(li_polyA_distant$gene, li_polyA_distant$SNP, sep="_")) #3 pairs
length(unique(str_split_fixed(eQTL_distant_same_SDP_gene, "_", 3)[,1])) # 1 genes regulated by same SDP in both tissue
eQTL_distant_specific_lv = setdiff(lv_polyA_distant$gene, li_polyA_distant$gene) # 31 genes
eQTL_distant_specific_li = setdiff(li_polyA_distant$gene, lv_polyA_distant$gene) # 11 genes
eQTL_distant_only_expressed_lv = intersect(lv_polyA_distant$gene, lv_specific) #13 genes with QTL only expressed in lv --> 0.68%
eQTL_distant_only_expressed_li = intersect(li_polyA_distant$gene, li_specific) #5 genes with QTL only expressed in li --> 0.70%


riboQTL_distant_same_gene = intersect(lv_ribo_distant$gene, li_ribo_distant$gene) #0 genes
riboQTL_distant_same_SDP_gene = intersect(paste(lv_ribo_distant$gene, lv_ribo_distant$SNP, sep="_"), paste(li_ribo_distant$gene, li_ribo_distant$SNP, sep="_")) #5 pairs
length(unique(str_split_fixed(riboQTL_distant_same_SDP_gene, "_", 3)[,1])) # 0 genes regulated by same SDP in both tissue
riboQTL_distant_specific_lv = setdiff(lv_ribo_distant$gene, li_ribo_distant$gene) # 45 genes
riboQTL_distant_specific_li = setdiff(li_ribo_distant$gene, lv_ribo_distant$gene) # 0 genes

riboQTL_distant_only_expressed_lv = intersect(lv_ribo_distant$gene, lv_specific) #8 genes with QTL only expressed in lv --> 0.42%
riboQTL_distant_only_expressed_li = intersect(li_ribo_distant$gene, li_specific) #0 genes with QTL only expressed in li --> 0%

teQTL_distant_same_gene = intersect(lv_te_distant$gene, li_te_distant$gene) #0 genes
teQTL_distant_same_SDP_gene = intersect(paste(lv_te_distant$gene, lv_te_distant$SNP, sep="_"), paste(li_te_distant$gene, li_te_distant$SNP, sep="_")) #0 pairs
length(unique(str_split_fixed(teQTL_distant_same_SDP_gene, "_", 3)[,1])) # 21 genes regulated by same SDP in both tissue

teQTL_distant_only_expressed_lv = intersect(lv_te_distant$gene, lv_specific) #16 genes with QTL only expressed in lv --> 0.84%
teQTL_distant_only_expressed_li = intersect(li_te_distant$gene, li_specific) #0 genes with QTL only expressed in li --> 0%

teQTL_distant_specific_lv = setdiff(lv_te_distant$gene, li_te_distant$gene) # 68 genes
teQTL_distant_specific_li = setdiff(li_te_distant$gene, lv_te_distant$gene) # 1 genes





#check for enrichment
# local eQTL
#enrichment of shared regulation?
fisher.test(matrix(c(126, 170, 434, 7892), nrow=2, byrow=T)) #2.2e-16
#enrichment of genes only expressed in lv
fisher.test(matrix(c(142, 560, 1767, 8062), nrow=2, byrow=T)) #not sign
#enrichment of genes only expressed in li
fisher.test(matrix(c(69, 296, 645, 8426), nrow=2, byrow=T)) #6.9e-13

# local riboQTL
#enrichment of shared regulation?
fisher.test(matrix(c(81, 240, 146, 7993), nrow=2, byrow=T)) #2.2e-16
#enrichment of genes only expressed in lv
fisher.test(matrix(c(82, 227, 1827, 8395), nrow=2, byrow=T)) #0.00022
#enrichment of genes only expressed in li
fisher.test(matrix(c(62, 321, 652, 8301), nrow=2, byrow=T)) #1.05e-08


# local teQTL
#enrichment of shared regulation?
fisher.test(matrix(c(21, 50, 39, 8512), nrow=2, byrow=T)) #2.2e-16
#enrichment of genes only expressed in lv
fisher.test(matrix(c(11, 60, 1898, 8562), nrow=2, byrow=T)) #not sign
#enrichment of genes only expressed in li
fisher.test(matrix(c(0, 71, 714, 8551), nrow=2, byrow=T)) #0.006



##real TE - not used for analysis!!!!

TE_polyA_RI_lv = read.table("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/TE/polyA_RI_lv_TE_deseqNorm_RNAANDRIBO.txt", sep="\t", header=T, stringsAsFactors=F)
TE_polyA_RI_li = read.table("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/TE/polyA_RI_li_TE_deseqNorm_RNAANDRIBO.txt", sep="\t", header=T, stringsAsFactors=F)

plot(density(log2(apply(TE_polyA_RI_lv[expressed_in_both,], 1, mean))), ylim=c(0,0.6))
lines(density(log2(apply(TE_polyA_RI_li[expressed_in_both,], 1, mean))), col="red")

##residuals - used for analysis
load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Input/Lv_Ribo_RNA_resid/Lv_Ribo_RNA_resid_input.RData")
resid_lv = pheno
load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Input/Liver_Ribo_RNA_resid/Liver_Ribo_RNA_resid_input.RData")
resid_li = pheno

plot(density(apply(resid_lv[expressed_in_both,], 1, function(x) mean(x, na.rm=T))), ylim=c(0,0.6))

lines(density(apply(resid_li[expressed_in_both,], 1, function(x) mean(x, na.rm=T))), col="red")





###







eQTL_distant_only_expressed_lv = intersect(lv_polyA_distant$gene, lv_specific) #13 genes with QTL only expressed in lv --> 0.68%
eQTL_distant_only_expressed_li = intersect(li_polyA_distant$gene, li_specific) #5 genes with QTL only expressed in li --> 0.70%

riboQTL_distant_only_expressed_lv = intersect(lv_ribo_distant$gene, lv_specific) #8 genes with QTL only expressed in lv --> 0.42%
riboQTL_distant_only_expressed_li = intersect(li_ribo_distant$gene, li_specific) #0 genes with QTL only expressed in li --> 0%

teQTL_distant_only_expressed_lv = intersect(lv_te_distant$gene, lv_specific) #16 genes with QTL only expressed in lv --> 0.84%
teQTL_distant_only_expressed_li = intersect(li_te_distant$gene, li_specific) #0 genes with QTL only expressed in li --> 0%

eQTL_local_only_expressed_lv = intersect(lv_polyA_local$gene, lv_specific) #142 genes with QTL only expressed in lv --> 7.44%
eQTL_local_only_expressed_li = intersect(li_polyA_local$gene, li_specific) #69 genes with QTL only expressed in li --> 8.26%

riboQTL_local_only_expressed_lv = intersect(lv_ribo_local$gene, lv_specific) #82 genes with QTL only expressed in lv --> 4.30%
riboQTL_local_only_expressed_li = intersect(li_ribo_local$gene, li_specific) #62 genes with QTL only expressed in li --> 8.68%

teQTL_local_only_expressed_lv = intersect(lv_te_local$gene, lv_specific) #11 genes with QTL only expressed in lv --> 0.58%
teQTL_local_only_expressed_li = intersect(li_te_local$gene, li_specific) #0 genes with QTL only expressed in li --> 0%






library(gplots)

lv_specific_associations = venn(list(polyA_local = unique(lv_polyA_local$gene), ribo_local = unique(lv_ribo_local$gene), te_local = unique(lv_te_local$gene), polyA_distant = unique(lv_polyA_distant$gene), ribo_distant = unique(lv_ribo_distant$gene), te_distant = unique(lv_te_distant$gene), all =  lv_specific), show.plot=F)


li_specific_associations = venn(list(polyA_local = unique(li_polyA_local$gene), ribo_local = unique(li_ribo_local$gene), te_local = unique(li_te_local$gene), polyA_distant = unique(li_polyA_distant$gene), ribo_distant = unique(li_ribo_distant$gene), te_distant = unique(li_te_distant$gene), all =  li_specific), show.plot=F)



