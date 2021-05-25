input_dir = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/permutations/"


lv_polyA_local = read.table(paste(input_dir, "Lv_polyA_esnps_cis.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)
lv_polyA_distant = read.table(paste(input_dir, "Lv_polyA_esnps_trans.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)


lv_ribo_local = read.table(paste(input_dir, "Lv_ribo_esnps_cis.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)
lv_ribo_distant = read.table(paste(input_dir, "Lv_ribo_esnps_trans.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)


lv_te_local = read.table(paste(input_dir, "Lv_RiboRNAresid_esnps_cis.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)
lv_te_distant = read.table(paste(input_dir, "Lv_RiboRNAresid_esnps_trans.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)


li_polyA_local = read.table(paste(input_dir, "Liver_polyA_esnps_cis.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)
li_polyA_distant = read.table(paste(input_dir, "Liver_polyA_esnps_trans.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)


li_ribo_local = read.table(paste(input_dir, "Liver_Ribo_esnps_cis.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)
li_ribo_distant = read.table(paste(input_dir, "Liver_Ribo_esnps_trans.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)


li_te_local = read.table(paste(input_dir, "Liver_RiboRNAresid_esnps_cis.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)
li_te_distant = read.table(paste(input_dir, "Liver_RiboRNAresid_esnps_trans.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)



#local polyA

lv_polyA_local$SNP_gene = paste(lv_polyA_local$SNP, lv_polyA_local$gene, sep="_")
li_polyA_local$SNP_gene = paste(li_polyA_local$SNP, li_polyA_local$gene, sep="_")

polyA_local_combined = data.frame(lv_polyA_local[, c("gene", "gene_name", "SNP_gene", "beta", "FDR_new", "empirical.p")], li_polyA_local[match(lv_polyA_local$SNP_gene, li_polyA_local$SNP_gene), c("beta", "FDR_new", "empirical.p")])
colnames(polyA_local_combined) = c("gene", "gene_name", "SNP_gene", "lv_eQTL_beta", "lv_eQTL_FDR_new", "lv_eQTL_emp.p","li_eQTL_beta", "li_eQTL_FDR_new", "li_eQTL_emp.p")
polyA_local_combined$ratio_beta = polyA_local_combined$lv_eQTL_beta / polyA_local_combined$li_eQTL_beta


#remove genes only tested in one tissue
polyA_local_combined.red = na.omit(polyA_local_combined)
polyA_local_combined.signInOne = polyA_local_combined[which((polyA_local_combined$lv_eQTL_FDR_new <= 0.1 & polyA_local_combined$lv_eQTL_emp.p<= 0.0015) | (polyA_local_combined$li_eQTL_FDR_new <= 0.1 & polyA_local_combined$li_eQTL_emp.p<= 0.0015)),]

##147 genes show tissue specific genetic regulation in lv or li
 
polyA_local_tissueSpecific = polyA_local_combined.signInOne[which(abs(polyA_local_combined.signInOne$ratio_beta) < 1/10 | abs(polyA_local_combined.signInOne$ratio_beta) > 10),]
length(unique(polyA_local_tissueSpecific[which(polyA_local_tissueSpecific$lv_eQTL_FDR_new <= 0.1 & polyA_local_tissueSpecific$lv_eQTL_emp.p<= 0.0015),]$gene_name))
#93 lv specific
length(unique(polyA_local_tissueSpecific[which(polyA_local_tissueSpecific$li_eQTL_FDR_new <= 0.1 & polyA_local_tissueSpecific$li_eQTL_emp.p<= 0.0015),]$gene_name))
#54 liver specific


load(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/filtered_genesets_ribotaper_fpkm.RData")
expressed_in_both = intersect(genesset_li, genesset_lv)

library(gProfileR)
set_base_url("http://biit.cs.ut.ee/gprofiler_archive2/r1477_e82_eg29/web")


local_eQTL_heart_go = gprofiler(query =  unique(polyA_local_tissueSpecific[which(polyA_local_tissueSpecific$lv_eQTL_FDR_new <= 0.1 & polyA_local_tissueSpecific$lv_eQTL_emp.p<= 0.0015),]$gene), custom_bg = expressed_in_both, organism = "rnorvegicus")
local_eQTL_liver_go = gprofiler(query =  unique(polyA_local_tissueSpecific[which(polyA_local_tissueSpecific$li_eQTL_FDR_new <= 0.1 & polyA_local_tissueSpecific$li_eQTL_emp.p<= 0.0015),]$gene), custom_bg = expressed_in_both, organism = "rnorvegicus")


polyA_local_combined.signInOne$tissue = "unspecific"
polyA_local_combined.signInOne[which((abs(polyA_local_combined.signInOne$ratio_beta) < 1/10 | abs(polyA_local_combined.signInOne$ratio_beta) > 10) & (polyA_local_combined.signInOne$lv_eQTL_FDR_new <= 0.1 & polyA_local_combined.signInOne$lv_eQTL_emp.p<= 0.0015)),]$tissue = "heart"
polyA_local_combined.signInOne[which((abs(polyA_local_combined.signInOne$ratio_beta) < 1/10 | abs(polyA_local_combined.signInOne$ratio_beta) > 10) & (polyA_local_combined.signInOne$li_eQTL_FDR_new <= 0.1 & polyA_local_combined.signInOne$li_eQTL_emp.p<= 0.0015)),]$tissue = "liver"


library(ggplot2)
pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/tissue_compare/test.pdf")
ggplot(data = polyA_local_combined.signInOne, aes(x=lv_eQTL_beta, y=li_eQTL_beta, col=tissue))+geom_point()+ylim(-3,3)+xlim(-3,3)+theme_classic()
hist(log2(abs(polyA_local_combined$ratio_beta)), breaks=60)
hist(log2(abs(polyA_local_combined[abs(polyA_local_combined$ratio_beta) <= 1/10,]$ratio_beta)), breaks=30, col="blue", add=T)
hist(log2(abs(polyA_local_combined[abs(polyA_local_combined$ratio_beta) >= 10,]$ratio_beta)), breaks=30, col="blue", add=T)
abline(v=log2(abs(10)))
abline(v=log2(abs(1/10)))
dev.off()

#local ribo

lv_ribo_local$SNP_gene = paste(lv_ribo_local$SNP, lv_ribo_local$gene, sep="_")
li_ribo_local$SNP_gene = paste(li_ribo_local$SNP, li_ribo_local$gene, sep="_")

ribo_local_combined = data.frame(lv_ribo_local[, c("gene", "gene_name", "SNP_gene", "beta", "FDR_new", "empirical.p")], li_ribo_local[match(lv_ribo_local$SNP_gene, li_ribo_local$SNP_gene), c("beta", "FDR_new", "empirical.p")])
colnames(ribo_local_combined) = c("gene", "gene_name", "SNP_gene", "lv_riboQTL_beta", "lv_riboQTL_FDR_new", "lv_riboQTL_emp.p","li_riboQTL_beta", "li_riboQTL_FDR_new", "li_riboQTL_emp.p")

ribo_local_combined$ratio_beta = ribo_local_combined$lv_riboQTL_beta / ribo_local_combined$li_riboQTL_beta
#remove genes only tested in one tissue
ribo_local_combined.red = na.omit(ribo_local_combined)
ribo_local_combined.signInOne = ribo_local_combined[which((ribo_local_combined$lv_riboQTL_FDR_new <= 0.1 & ribo_local_combined$lv_riboQTL_emp.p<= 0.0015) | (ribo_local_combined$li_riboQTL_FDR_new <= 0.1 & ribo_local_combined$li_riboQTL_emp.p<= 0.0015)),]

##149 genes show tissue specific genetic regulation in lv or li
ribo_local_tissueSpecific = ribo_local_combined.signInOne[which(abs(ribo_local_combined.signInOne$ratio_beta) < 1/10 | abs(ribo_local_combined.signInOne$ratio_beta) > 10),]
length(unique(ribo_local_tissueSpecific[which(ribo_local_tissueSpecific$lv_riboQTL_FDR_new <= 0.1 & ribo_local_tissueSpecific$lv_riboQTL_emp.p<= 0.0015),]$gene_name))
#38 lv specific
length(unique(ribo_local_tissueSpecific[which(ribo_local_tissueSpecific$li_riboQTL_FDR_new <= 0.1 & ribo_local_tissueSpecific$li_riboQTL_emp.p<= 0.0015),]$gene_name))
#66 liver specific




ribo_local_combined.signInOne$tissue = "unspecific"
ribo_local_combined.signInOne[which((abs(ribo_local_combined.signInOne$ratio_beta) < 1/5 | abs(ribo_local_combined.signInOne$ratio_beta) > 5) & (ribo_local_combined.signInOne$lv_riboQTL_FDR_new <= 0.1 & ribo_local_combined.signInOne$lv_riboQTL_emp.p<= 0.0015)),]$tissue = "heart"
ribo_local_combined.signInOne[which((abs(ribo_local_combined.signInOne$ratio_beta) < 1/5 | abs(ribo_local_combined.signInOne$ratio_beta) > 5) & (ribo_local_combined.signInOne$li_riboQTL_FDR_new <= 0.1 & ribo_local_combined.signInOne$li_riboQTL_emp.p<= 0.0015)),]$tissue = "liver"


library(ggplot2)
pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/tissue_compare/test_ribo.pdf")
ggplot(data = ribo_local_combined.signInOne, aes(x=lv_riboQTL_beta, y=li_riboQTL_beta, col=tissue))+geom_point()+ylim(-3,3)+xlim(-3,3)+theme_classic()
hist(log2(abs(ribo_local_combined$ratio_beta)), breaks=60)
hist(log2(abs(ribo_local_combined[abs(ribo_local_combined$ratio_beta) <= 1/10,]$ratio_beta)), breaks=30, col="blue", add=T)
hist(log2(abs(ribo_local_combined[abs(ribo_local_combined$ratio_beta) >= 10,]$ratio_beta)), breaks=30, col="blue", add=T)
abline(v=log2(abs(10)))
abline(v=log2(abs(1/10)))
dev.off()







#distant polyA
lv_polyA_distant$SNP_gene = paste(lv_polyA_distant$SNP, lv_polyA_distant$gene, sep="_")
li_polyA_distant$SNP_gene = paste(li_polyA_distant$SNP, li_polyA_distant$gene, sep="_")


polyA_distant_combined = data.frame(lv_polyA_distant[, c("gene_name", "SNP_gene", "beta", "FDR_new", "empirical.p")], li_polyA_distant[match(lv_polyA_distant$SNP_gene, li_polyA_distant$SNP_gene), c("beta", "FDR_new", "empirical.p")])
colnames(polyA_distant_combined) = c("gene_name", "SNP_gene", "lv_eQTL_beta", "lv_eQTL_FDR_new", "lv_eQTL_emp.p","li_eQTL_beta", "li_eQTL_FDR_new", "li_eQTL_emp.p")

polyA_distant_combined$ratio_beta = polyA_distant_combined$lv_eQTL_beta / polyA_distant_combined$li_eQTL_beta
#remove genes only tested in one tissue
polyA_distant_combined.red = na.omit(polyA_distant_combined)
polyA_distant_combined.signInOne = polyA_distant_combined[which((polyA_distant_combined$lv_eQTL_FDR_new <= 0.2 & polyA_distant_combined$lv_eQTL_emp.p<= 0.0015) | (polyA_distant_combined$li_eQTL_FDR_new <= 0.2 & polyA_distant_combined$li_eQTL_emp.p<= 0.0015)),]

##11 genes show tissue specific genetic regulation in lv or li
polyA_distant_tissueSpecific = polyA_distant_combined.signInOne[which(abs(polyA_distant_combined.signInOne$ratio_beta) < 1/10 | abs(polyA_distant_combined.signInOne$ratio_beta) > 10),]
polyA_distant_tissueSpecific[which(polyA_distant_tissueSpecific$lv_eQTL_FDR_new <= 0.2 & polyA_distant_tissueSpecific$lv_eQTL_emp.p<= 0.0015),]
#5 lv specific ("Il2rb"    "Hist1h1c" "Myo5b"    "Dph5"     "Mbp"      "Pten"     "Copa"     "Ppm1f"    "Taco1")
polyA_distant_tissueSpecific[which(polyA_distant_tissueSpecific$li_eQTL_FDR_new <= 0.2 & polyA_distant_tissueSpecific$li_eQTL_emp.p<= 0.0015),]
#3 liver specific (gas7, S100a10)




polyA_distant_combined.signInOne$tissue = "unspecific"
polyA_distant_combined.signInOne[which((abs(polyA_distant_combined.signInOne$ratio_beta) < 1/5 | abs(polyA_distant_combined.signInOne$ratio_beta) > 5) & (polyA_distant_combined.signInOne$lv_eQTL_FDR_new <= 0.2 & polyA_distant_combined.signInOne$lv_eQTL_emp.p<= 0.0015)),]$tissue = "heart"
polyA_distant_combined.signInOne[which((abs(polyA_distant_combined.signInOne$ratio_beta) < 1/5 | abs(polyA_distant_combined.signInOne$ratio_beta) > 5) & (polyA_distant_combined.signInOne$li_eQTL_FDR_new <= 0.2 & polyA_distant_combined.signInOne$li_eQTL_emp.p<= 0.0015)),]$tissue = "liver"


library(ggplot2)
pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/tissue_compare/test_distant.pdf")
ggplot(data = polyA_distant_combined.signInOne, aes(x=lv_eQTL_beta, y=li_eQTL_beta, col=tissue))+geom_point()+ylim(-3,3)+xlim(-3,3)
dev.off()




#distant ribo

lv_ribo_distant$SNP_gene = paste(lv_ribo_distant$SNP, lv_ribo_distant$gene, sep="_")
li_ribo_distant$SNP_gene = paste(li_ribo_distant$SNP, li_ribo_distant$gene, sep="_")

ribo_distant_combined = data.frame(lv_ribo_distant[, c("gene_name", "SNP_gene", "beta", "FDR_new", "empirical.p")], li_ribo_distant[match(lv_ribo_distant$SNP_gene, li_ribo_distant$SNP_gene), c("beta", "FDR_new", "empirical.p")])
colnames(ribo_distant_combined) = c("gene_name", "SNP_gene", "lv_riboQTL_beta", "lv_riboQTL_FDR_new", "lv_riboQTL_emp.p","li_riboQTL_beta", "li_riboQTL_FDR_new", "li_riboQTL_emp.p")

ribo_distant_combined$ratio_beta = ribo_distant_combined$lv_riboQTL_beta / ribo_distant_combined$li_riboQTL_beta
#remove genes only tested in one tissue
ribo_distant_combined.red = na.omit(ribo_distant_combined)
ribo_distant_combined.signInOne = ribo_distant_combined[which((ribo_distant_combined$lv_riboQTL_FDR_new <= 0.2 & ribo_distant_combined$lv_riboQTL_emp.p<= 0.0015) | (ribo_distant_combined$li_riboQTL_FDR_new <= 0.2 & ribo_distant_combined$li_riboQTL_emp.p<= 0.0015)),]

##149 genes show tissue specific genetic regulation in lv or li
ribo_distant_tissueSpecific = ribo_distant_combined.signInOne[which(abs(ribo_distant_combined.signInOne$ratio_beta) < 1/10 | abs(ribo_distant_combined.signInOne$ratio_beta) > 10),]
ribo_distant_tissueSpecific[which(ribo_distant_tissueSpecific$lv_riboQTL_FDR_new <= 0.2 & ribo_distant_tissueSpecific$lv_riboQTL_emp.p<= 0.0015),]
#16 lv specific
ribo_distant_tissueSpecific[which(ribo_distant_tissueSpecific$li_riboQTL_FDR_new <= 0.2 & ribo_distant_tissueSpecific$li_riboQTL_emp.p<= 0.0015),]
#0 liver specific

save.image("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/tissue_compare/tissue_specific_QTLs.RData")


















