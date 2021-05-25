res_dir="/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/"


load(file = paste(res_dir, "Result_tables/20180213_Lv_polyA/eqtls.RData", sep=""))
qtls_lv_rna = actual.eqtls
load(file=paste(res_dir, "Result_tables/20180213_Lv_Ribo/eqtls.RData", sep=""))
qtls_lv_ribo = actual.eqtls
load(file=paste(res_dir, "Result_tables/20180213_Lv_Ribo_RNA_resid/eqtls.RData", sep=""))
qtls_lv_te = actual.eqtls

load(file = paste(res_dir, "Result_tables/20201117_Lv_polyA_sva/eqtls.RData", sep=""))
qtls_lv_rna_sva = actual.eqtls
load(file=paste(res_dir, "Result_tables/20201117_Lv_Ribo_sva/eqtls.RData", sep=""))
qtls_lv_ribo_sva = actual.eqtls
load(file=paste(res_dir, "Result_tables/20201117_Lv_Ribo_RNA_resid_sva/eqtls.RData", sep=""))
qtls_lv_te_sva = actual.eqtls


qtls_lv_rna = rbind(qtls_lv_rna[["cis"]], qtls_lv_rna[["trans"]])
qtls_lv_rna_sva = rbind(qtls_lv_rna_sva[["cis"]], qtls_lv_rna_sva[["trans"]])
qtls_lv_ribo = rbind(qtls_lv_ribo[["cis"]], qtls_lv_ribo[["trans"]])
qtls_lv_ribo_sva = rbind(qtls_lv_ribo_sva[["cis"]], qtls_lv_ribo_sva[["trans"]])
qtls_lv_te = rbind(qtls_lv_te[["cis"]], qtls_lv_te[["trans"]])
qtls_lv_te_sva = rbind(qtls_lv_te_sva[["cis"]], qtls_lv_te_sva[["trans"]])



qtls_lv_rna$key = paste(qtls_lv_rna$SNP, qtls_lv_rna$gene, sep="_")
qtls_lv_rna_sva$key = paste(qtls_lv_rna_sva$SNP, qtls_lv_rna_sva$gene, sep="_")

RNA_data = data.frame(key1 = qtls_lv_rna$key, SNP = qtls_lv_rna$SNP, gene = qtls_lv_rna$gene, beta_old = qtls_lv_rna$beta, 
                      pvalue_old = qtls_lv_rna$`p-value`, beta_sva = qtls_lv_rna_sva[match(qtls_lv_rna$key, qtls_lv_rna_sva$key), "beta"], 
                      pvalue_sva = qtls_lv_rna_sva[match(qtls_lv_rna$key, qtls_lv_rna_sva$key), "p-value"],
                      key2 = qtls_lv_rna_sva[match(qtls_lv_rna$key, qtls_lv_rna_sva$key), "key"] )

lv_rna_local = read.table("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/local_eQTLs_Lv_polyA_sign.txt", sep="\t", header=T, stringsAsFactors = F)
lv_rna_distant = read.table("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/distant_eQTLs_Lv_polyA_sign.txt", sep="\t", header=T, stringsAsFactors = F)
lv_rna_local$key = paste(lv_rna_local$SNP, lv_rna_local$gene, sep="_")
lv_rna_distant$key = paste(lv_rna_distant$SNP, lv_rna_distant$gene, sep="_")

RNA_data_filtLocal = RNA_data[match(lv_rna_local$key, RNA_data$key1),]
RNA_data_filtDistant = RNA_data[match(lv_rna_distant$key, RNA_data$key1),]

cor.test(RNA_data_filtLocal$beta, RNA_data_filtLocal$beta_old) #0.996
cor.test(RNA_data_filtLocal$p.value, RNA_data_filtLocal$pvalue_old, method = "spearman") #0.84

cor.test(RNA_data_filtDistant$beta, RNA_data_filtDistant$beta_old) # 0.999
cor.test(RNA_data_filtDistant$p.value, RNA_data_filtDistant$pvalue_old,method = "spearman") #0.84


qtls_lv_ribo$key = paste(qtls_lv_ribo$SNP, qtls_lv_ribo$gene, sep="_")
qtls_lv_ribo_sva$key = paste(qtls_lv_ribo_sva$SNP, qtls_lv_ribo_sva$gene, sep="_")

Ribo_data = data.frame(key1 = qtls_lv_ribo$key, SNP = qtls_lv_ribo$SNP, gene_name = qtls_lv_ribo$gene_name, gene = qtls_lv_ribo$gene, beta_old = qtls_lv_ribo$beta, 
                      pvalue_old = qtls_lv_ribo$`p-value`, beta_sva = qtls_lv_ribo_sva[match(qtls_lv_ribo$key, qtls_lv_ribo_sva$key), "beta"], 
                      pvalue_sva = qtls_lv_ribo_sva[match(qtls_lv_ribo$key, qtls_lv_ribo_sva$key), "p-value"],
                      key2 = qtls_lv_ribo_sva[match(qtls_lv_ribo$key, qtls_lv_ribo_sva$key), "key"] )

lv_ribo_local = read.table("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/local_eQTLs_Lv_Ribo_sign.txt", sep="\t", header=T, stringsAsFactors = F)
lv_ribo_distant = read.table("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/distant_eQTLs_Lv_Ribo_sign.txt", sep="\t", header=T, stringsAsFactors = F)
lv_ribo_local$key = paste(lv_ribo_local$SNP, lv_ribo_local$gene, sep="_")
lv_ribo_distant$key = paste(lv_ribo_distant$SNP, lv_ribo_distant$gene, sep="_")

Ribo_data_filtLocal = Ribo_data[match(lv_ribo_local$key, Ribo_data$key1),]
Ribo_data_filtDistant = Ribo_data[match(lv_ribo_distant$key, Ribo_data$key1),]

cor.test(Ribo_data_filtLocal$beta, Ribo_data_filtLocal$beta_old) #0.9777
cor.test(Ribo_data_filtLocal$p.value, Ribo_data_filtLocal$pvalue_old, method = "spearman") #0.544

cor.test(Ribo_data_filtDistant$beta, Ribo_data_filtDistant$beta_old) # 0.989
cor.test(Ribo_data_filtDistant$p.value, Ribo_data_filtDistant$pvalue_old,method = "spearman") #0.886
write.table(Ribo_data_filtDistant, file = "Desktop/Ribo_data_filtDistant_SVA.txt", sep="\t", col.names=T, row.names = F, quote = F)
save.image("Desktop/MatrixEQTL_withSVA.RData")





qtls_lv_te$key = paste(qtls_lv_te$SNP, qtls_lv_te$gene, sep="_")
qtls_lv_te_sva$key = paste(qtls_lv_te_sva$SNP, qtls_lv_te_sva$gene, sep="_")

te_data = data.frame(key1 = qtls_lv_te$key, SNP = qtls_lv_te$SNP, gene_name = qtls_lv_te$gene_name, gene = qtls_lv_te$gene, beta_old = qtls_lv_te$beta, 
                       pvalue_old = qtls_lv_te$`p-value`, beta_sva = qtls_lv_te_sva[match(qtls_lv_te$key, qtls_lv_te_sva$key), "beta"], 
                       pvalue_sva = qtls_lv_te_sva[match(qtls_lv_te$key, qtls_lv_te_sva$key), "p-value"],
                       key2 = qtls_lv_te_sva[match(qtls_lv_te$key, qtls_lv_te_sva$key), "key"] )

lv_te_local = read.table("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/local_eQTLs_Lv_RiboRNAresid_sign.txt", sep="\t", header=T, stringsAsFactors = F)
lv_te_distant = read.table("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/distant_eQTLs_Lv_RiboRNAresid_sign.txt", sep="\t", header=T, stringsAsFactors = F)
lv_te_local$key = paste(lv_te_local$SNP, lv_te_local$gene, sep="_")
lv_te_distant$key = paste(lv_te_distant$SNP, lv_te_distant$gene, sep="_")

te_data_filtLocal = te_data[match(lv_te_local$key, te_data$key1),]
te_data_filtDistant = te_data[match(lv_te_distant$key, te_data$key1),]

cor.test(te_data_filtLocal$beta, te_data_filtLocal$beta_old) #0.983397
cor.test(te_data_filtLocal$p.value, te_data_filtLocal$pvalue_old, method = "spearman") #0.544

cor.test(te_data_filtDistant$beta, te_data_filtDistant$beta_old) # 0.989
cor.test(te_data_filtDistant$p.value, te_data_filtDistant$pvalue_old,method = "spearman") #0.524
write.table(te_data_filtDistant, file = "Desktop/TE_data_filtDistant_SVA.txt", sep="\t", col.names=T, row.names = F, quote = F)
save.image("Desktop/MatrixEQTL_withSVA.RData")


library(ggplot2)
library(gridExtra)
pdf("Desktop/Scatter_QTLMapping_CompareSVA_noSVA.pdf", width = 15, height=10)
a=ggplot(data = RNA_data_filtDistant, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant eQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.5))+
  scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.5))+geom_text(aes(-1.3, 1.5 , label=paste0("R = ", round(cor.test(RNA_data_filtDistant$beta, RNA_data_filtDistant$beta_old)$estimate, 3)), vjust=0))
b = ggplot(data = Ribo_data_filtDistant, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant riboQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.7))+
  scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.7))+geom_text(aes(-1.3, 1.7 , label=paste0("R = ", round(cor.test(Ribo_data_filtDistant$beta, Ribo_data_filtDistant$beta_old)$estimate, 3)), vjust=0))
c=ggplot(data = te_data_filtDistant, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant teQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.5))+
  scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.5))+geom_text(aes(-1.3, 1.5 , label=paste0("R = ", round(cor.test(te_data_filtDistant$beta, te_data_filtDistant$beta_old)$estimate, 3)), vjust=0))

d=ggplot(data = RNA_data_filtLocal, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local eQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-4,4,1), limits = c(-4,4))+
  scale_y_continuous(breaks = seq(-4,4,1), limits = c(-4,4)) +geom_text(aes(-3.5, 4 , label=paste0("R = ", round(cor.test(RNA_data_filtLocal$beta, RNA_data_filtLocal$beta_old)$estimate, 3)), vjust=0))

e=ggplot(data = Ribo_data_filtLocal, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local riboQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-3,5.5,1), limits = c(-3,5))+
  scale_y_continuous(breaks =  seq(-3,5,1), limits = c(-3,5)) +geom_text(aes(-2.5, 5 , label=paste0("R = ", round(cor.test(Ribo_data_filtLocal$beta, Ribo_data_filtLocal$beta_old)$estimate, 3)), vjust=0))

f=ggplot(data = te_data_filtLocal, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local teQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-1.5,1.5,0.5), limits=c(-2,1.6))+
  scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits=c(-2,1.6)) +geom_text(aes(-1.7, 1.6 , label=paste0("R = ", round(cor.test(te_data_filtLocal$beta, te_data_filtLocal$beta_old)$estimate, 3)), vjust=0))
grid.arrange(d,e,f,a,b,c, ncol=3, nrow=2)
dev.off()

pdf("Desktop/Scatter_QTLMapping_CompareSVA_noSVA_pvalues.pdf", width = 15, height=10)
a=ggplot(data = RNA_data_filtDistant, aes(x = log10(p.value), y = log10(pvalue_old))) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant eQTLs Left Ventricle") +
  xlab("log10(p-values) using relatedness as covariate")+ylab("log10(pvalues) without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-15,1,2), limits = c(-15, 1))+
  scale_y_continuous(breaks = seq(-15,0,2), limits = c(-15, 0))+geom_text(aes(-14, 0 , label=paste0("R = ", round(cor.test(log10(RNA_data_filtDistant$p.value), log10(RNA_data_filtDistant$pvalue_old))$estimate, 3)), vjust=0))

  
b = ggplot(data = Ribo_data_filtDistant, aes(x = log10(p.value), y = log10(pvalue_old))) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant riboQTLs Left Ventricle") +
  xlab("log10(p-values) using relatedness as covariate")+ylab("log10(pvalues) without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-16,1,2), limits = c(-16, 1))+
  scale_y_continuous(breaks = seq(-16,0,2), limits = c(-16, 0))+geom_text(aes(-15, 0 , label=paste0("R = ", round(cor.test(log10(Ribo_data_filtDistant$p.value), log10(Ribo_data_filtDistant$pvalue_old))$estimate, 3)), vjust=0))

  
c=ggplot(data = te_data_filtDistant, aes(x = log10(p.value), y = log10(pvalue_old))) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant teQTLs Left Ventricle") +
  xlab("log10(p-values) using relatedness as covariate")+ylab("log10(pvalues) without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-18,1,2), limits = c(-18, 1))+
  scale_y_continuous(breaks = seq(-18,0,2), limits = c(-18, 0))+geom_text(aes(-16.5, 0 , label=paste0("R = ", round(cor.test(log10(te_data_filtDistant$p.value), log10(te_data_filtDistant$pvalue_old))$estimate, 3)), vjust=0))

d=ggplot(data = RNA_data_filtLocal, aes(x = log10(p.value), y = log10(pvalue_old))) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local eQTLs Left Ventricle") +
  xlab("log10(p-values) using relatedness as covariate")+ylab("log10(pvalues) without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-30,1,5), limits = c(-30, 1))+
  scale_y_continuous(breaks = seq(-30,0,5), limits = c(-30, 0))+geom_text(aes(-28, 0 , label=paste0("R = ", round(cor.test(log10(RNA_data_filtLocal$p.value), log10(RNA_data_filtLocal$pvalue_old))$estimate, 3)), vjust=0))

e=ggplot(data = Ribo_data_filtLocal, aes(x = log10(p.value), y = log10(pvalue_old))) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local riboQTLs Left Ventricle") +
  xlab("log10(p-values) using relatedness as covariate")+ylab("log10(pvalues) without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-30,1,5), limits = c(-30, 1))+
  scale_y_continuous(breaks = seq(-30,0,5), limits = c(-30, 0))+geom_text(aes(-28, 0 , label=paste0("R = ", round(cor.test(log10(Ribo_data_filtLocal$p.value), log10(Ribo_data_filtLocal$pvalue_old))$estimate, 3)), vjust=0))

f=ggplot(data = te_data_filtLocal, aes(x = log10(p.value), y = log10(pvalue_old))) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local teQTLs Left Ventricle") +
  xlab("log10(p-values) using relatedness as covariate")+ylab("log10(pvalues) without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-20,1,2.5), limits = c(-20, 1))+
  scale_y_continuous(breaks = seq(-20,0,2.5), limits = c(-20, 0))+geom_text(aes(-19, 0 , label=paste0("R = ", round(cor.test(log10(te_data_filtLocal$p.value), log10(te_data_filtLocal$pvalue_old))$estimate, 3)), vjust=0))

  
  grid.arrange(d,e,f,a,b,c, ncol=3, nrow=2)
dev.off()

###liver

load(file = paste(res_dir, "Result_tables/20180213_Liver_polyA/eqtls.RData", sep=""))
qtls_lv_rna = actual.eqtls
load(file=paste(res_dir, "Result_tables/20180213_Liver_Ribo/eqtls.RData", sep=""))
qtls_lv_ribo = actual.eqtls
load(file=paste(res_dir, "Result_tables/20180213_Liver_Ribo_RNA_resid/eqtls.RData", sep=""))
qtls_lv_te = actual.eqtls

load(file = paste(res_dir, "Result_tables/20201119_Liver_polyA_sva/eqtls.RData", sep=""))
qtls_liver_rna_sva = actual.eqtls
load(file=paste(res_dir, "Result_tables/20201119_Liver_Ribo_sva/eqtls.RData", sep=""))
qtls_liver_ribo_sva = actual.eqtls
load(file=paste(res_dir, "Result_tables/20201119_Liver_Ribo_RNA_resid_sva/eqtls.RData", sep=""))
qtls_liver_te_sva = actual.eqtls


qtls_liver_rna = rbind(qtls_liver_rna[["cis"]], qtls_liver_rna[["trans"]])
qtls_liver_rna_sva = rbind(qtls_liver_rna_sva[["cis"]], qtls_liver_rna_sva[["trans"]])
qtls_liver_ribo = rbind(qtls_liver_ribo[["cis"]], qtls_liver_ribo[["trans"]])
qtls_liver_ribo_sva = rbind(qtls_liver_ribo_sva[["cis"]], qtls_liver_ribo_sva[["trans"]])
qtls_liver_te = rbind(qtls_liver_te[["cis"]], qtls_liver_te[["trans"]])
qtls_liver_te_sva = rbind(qtls_liver_te_sva[["cis"]], qtls_liver_te_sva[["trans"]])



qtls_liver_rna$key = paste(qtls_liver_rna$SNP, qtls_liver_rna$gene, sep="_")
qtls_liver_rna_sva$key = paste(qtls_liver_rna_sva$SNP, qtls_liver_rna_sva$gene, sep="_")

RNA_data = data.frame(key1 = qtls_liver_rna$key, SNP = qtls_liver_rna$SNP, gene = qtls_liver_rna$gene, beta_old = qtls_liver_rna$beta, 
                      pvalue_old = qtls_liver_rna$`p-value`, beta_sva = qtls_liver_rna_sva[match(qtls_liver_rna$key, qtls_liver_rna_sva$key), "beta"], 
                      pvalue_sva = qtls_liver_rna_sva[match(qtls_liver_rna$key, qtls_liver_rna_sva$key), "p-value"],
                      key2 = qtls_liver_rna_sva[match(qtls_liver_rna$key, qtls_liver_rna_sva$key), "key"] )

liver_rna_local = read.table("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/local_eQTLs_Liver_polyA_sign.txt", sep="\t", header=T, stringsAsFactors = F)
liver_rna_distant = read.table("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/distant_eQTLs_Liver_polyA_sign.txt", sep="\t", header=T, stringsAsFactors = F)
liver_rna_local$key = paste(liver_rna_local$SNP, liver_rna_local$gene, sep="_")
liver_rna_distant$key = paste(liver_rna_distant$SNP, liver_rna_distant$gene, sep="_")

RNA_data_filtLocal = RNA_data[match(liver_rna_local$key, RNA_data$key1),]
RNA_data_filtDistant = RNA_data[match(liver_rna_distant$key, RNA_data$key1),]

cor.test(RNA_data_filtLocal$beta, RNA_data_filtLocal$beta_old) #0.996
cor.test(RNA_data_filtLocal$p.value, RNA_data_filtLocal$pvalue_old, method = "spearman") #0.84

cor.test(RNA_data_filtDistant$beta, RNA_data_filtDistant$beta_old) # 0.999
cor.test(RNA_data_filtDistant$p.value, RNA_data_filtDistant$pvalue_old,method = "spearman") #0.84


qtls_liver_ribo$key = paste(qtls_liver_ribo$SNP, qtls_liver_ribo$gene, sep="_")
qtls_liver_ribo_sva$key = paste(qtls_liver_ribo_sva$SNP, qtls_liver_ribo_sva$gene, sep="_")

Ribo_data = data.frame(key1 = qtls_liver_ribo$key, SNP = qtls_liver_ribo$SNP, gene_name = qtls_liver_ribo$gene_name, gene = qtls_liver_ribo$gene, beta_old = qtls_liver_ribo$beta, 
                       pvalue_old = qtls_liver_ribo$`p-value`, beta_sva = qtls_liver_ribo_sva[match(qtls_liver_ribo$key, qtls_liver_ribo_sva$key), "beta"], 
                       pvalue_sva = qtls_liver_ribo_sva[match(qtls_liver_ribo$key, qtls_liver_ribo_sva$key), "p-value"],
                       key2 = qtls_liver_ribo_sva[match(qtls_liver_ribo$key, qtls_liver_ribo_sva$key), "key"] )

liver_ribo_local = read.table("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/local_eQTLs_Liver_Ribo_sign.txt", sep="\t", header=T, stringsAsFactors = F)
liver_ribo_distant = read.table("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/distant_eQTLs_Liver_Ribo_sign.txt", sep="\t", header=T, stringsAsFactors = F)
liver_ribo_local$key = paste(liver_ribo_local$SNP, liver_ribo_local$gene, sep="_")
liver_ribo_distant$key = paste(liver_ribo_distant$SNP, liver_ribo_distant$gene, sep="_")

Ribo_data_filtLocal = Ribo_data[match(liver_ribo_local$key, Ribo_data$key1),]
Ribo_data_filtDistant = Ribo_data[match(liver_ribo_distant$key, Ribo_data$key1),]

cor.test(Ribo_data_filtLocal$beta, Ribo_data_filtLocal$beta_old) #0.9777
cor.test(Ribo_data_filtLocal$p.value, Ribo_data_filtLocal$pvalue_old, method = "spearman") #0.544

cor.test(Ribo_data_filtDistant$beta, Ribo_data_filtDistant$beta_old) # 0.989
cor.test(Ribo_data_filtDistant$p.value, Ribo_data_filtDistant$pvalue_old,method = "spearman") #0.886
write.table(Ribo_data_filtDistant, file = "Desktop/Ribo_data_filtDistant_SVA_liver.txt", sep="\t", col.names=T, row.names = F, quote = F)
save.image("Desktop/MatrixEQTL_withSVA_liver.RData")





qtls_liver_te$key = paste(qtls_liver_te$SNP, qtls_liver_te$gene, sep="_")
qtls_liver_te_sva$key = paste(qtls_liver_te_sva$SNP, qtls_liver_te_sva$gene, sep="_")

te_data = data.frame(key1 = qtls_liver_te$key, SNP = qtls_liver_te$SNP, gene_name = qtls_liver_te$gene_name, gene = qtls_liver_te$gene, beta_old = qtls_liver_te$beta, 
                     pvalue_old = qtls_liver_te$`p-value`, beta_sva = qtls_liver_te_sva[match(qtls_liver_te$key, qtls_liver_te_sva$key), "beta"], 
                     pvalue_sva = qtls_liver_te_sva[match(qtls_liver_te$key, qtls_liver_te_sva$key), "p-value"],
                     key2 = qtls_liver_te_sva[match(qtls_liver_te$key, qtls_liver_te_sva$key), "key"] )

liver_te_local = read.table("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/local_eQTLs_Liver_RiboRNAresid_sign.txt", sep="\t", header=T, stringsAsFactors = F)
liver_te_distant = read.table("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/distant_eQTLs_Liver_RiboRNAresid_sign.txt", sep="\t", header=T, stringsAsFactors = F)
liver_te_local$key = paste(liver_te_local$SNP, liver_te_local$gene, sep="_")
liver_te_distant$key = paste(liver_te_distant$SNP, liver_te_distant$gene, sep="_")

te_data_filtLocal = te_data[match(liver_te_local$key, te_data$key1),]
te_data_filtDistant = te_data[match(liver_te_distant$key, te_data$key1),]

cor.test(te_data_filtLocal$beta, te_data_filtLocal$beta_old) #0.983397
cor.test(te_data_filtLocal$p.value, te_data_filtLocal$pvalue_old, method = "spearman") #0.544

cor.test(te_data_filtDistant$beta, te_data_filtDistant$beta_old) # 0.989
cor.test(te_data_filtDistant$p.value, te_data_filtDistant$pvalue_old,method = "spearman") #0.524
write.table(te_data_filtDistant, file = "Desktop/TE_data_filtDistant_SVA_liver.txt", sep="\t", col.names=T, row.names = F, quote = F)
save.image("Desktop/MatrixEQTL_withSVA_liver.RData")


library(ggplot2)
library(gridExtra)
pdf("Desktop/Scatter_QTLMapping_CompareSVA_noSVA_liver.pdf", width = 15, height=10)
a=ggplot(data = RNA_data_filtDistant, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant eQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.5))+
  scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.5))+geom_text(aes(-1.3, 1.5 , label=paste0("R = ", round(cor.test(RNA_data_filtDistant$beta, RNA_data_filtDistant$beta_old)$estimate, 3)), vjust=0))
b = ggplot(data = Ribo_data_filtDistant, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant riboQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.7))+
  scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.7))+geom_text(aes(-1.3, 1.7 , label=paste0("R = ", round(cor.test(Ribo_data_filtDistant$beta, Ribo_data_filtDistant$beta_old)$estimate, 3)), vjust=0))
c=ggplot(data = te_data_filtDistant, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant teQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.5))+
  scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.5))+geom_text(aes(-1.3, 1.5 , label=paste0("R = ", round(cor.test(te_data_filtDistant$beta, te_data_filtDistant$beta_old)$estimate, 3)), vjust=0))

d=ggplot(data = RNA_data_filtLocal, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local eQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-4,4,1), limits = c(-4,4))+
  scale_y_continuous(breaks = seq(-4,4,1), limits = c(-4,4)) +geom_text(aes(-3.5, 4 , label=paste0("R = ", round(cor.test(RNA_data_filtLocal$beta, RNA_data_filtLocal$beta_old)$estimate, 3)), vjust=0))

e=ggplot(data = Ribo_data_filtLocal, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local riboQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-3,5.5,1), limits = c(-3,5))+
  scale_y_continuous(breaks =  seq(-3,5,1), limits = c(-3,5)) +geom_text(aes(-2.5, 5 , label=paste0("R = ", round(cor.test(Ribo_data_filtLocal$beta, Ribo_data_filtLocal$beta_old)$estimate, 3)), vjust=0))

f=ggplot(data = te_data_filtLocal, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local teQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-1.5,1.5,0.5), limits=c(-2,1.6))+
  scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits=c(-2,1.6)) +geom_text(aes(-1.7, 1.6 , label=paste0("R = ", round(cor.test(te_data_filtLocal$beta, te_data_filtLocal$beta_old)$estimate, 3)), vjust=0))
grid.arrange(d,e,f,a,b,c, ncol=3, nrow=2)
dev.off()

pdf("Desktop/Scatter_QTLMapping_CompareSVA_noSVA_pvalues.pdf", width = 15, height=10)
a=ggplot(data = RNA_data_filtDistant, aes(x = log10(p.value), y = log10(pvalue_old))) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant eQTLs Left Ventricle") +
  xlab("log10(p-values) using relatedness as covariate")+ylab("log10(pvalues) without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-15,1,2), limits = c(-15, 1))+
  scale_y_continuous(breaks = seq(-15,0,2), limits = c(-15, 0))+geom_text(aes(-14, 0 , label=paste0("R = ", round(cor.test(log10(RNA_data_filtDistant$p.value), log10(RNA_data_filtDistant$pvalue_old))$estimate, 3)), vjust=0))


b = ggplot(data = Ribo_data_filtDistant, aes(x = log10(p.value), y = log10(pvalue_old))) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant riboQTLs Left Ventricle") +
  xlab("log10(p-values) using relatedness as covariate")+ylab("log10(pvalues) without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-16,1,2), limits = c(-16, 1))+
  scale_y_continuous(breaks = seq(-16,0,2), limits = c(-16, 0))+geom_text(aes(-15, 0 , label=paste0("R = ", round(cor.test(log10(Ribo_data_filtDistant$p.value), log10(Ribo_data_filtDistant$pvalue_old))$estimate, 3)), vjust=0))


c=ggplot(data = te_data_filtDistant, aes(x = log10(p.value), y = log10(pvalue_old))) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant teQTLs Left Ventricle") +
  xlab("log10(p-values) using relatedness as covariate")+ylab("log10(pvalues) without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-18,1,2), limits = c(-18, 1))+
  scale_y_continuous(breaks = seq(-18,0,2), limits = c(-18, 0))+geom_text(aes(-16.5, 0 , label=paste0("R = ", round(cor.test(log10(te_data_filtDistant$p.value), log10(te_data_filtDistant$pvalue_old))$estimate, 3)), vjust=0))

d=ggplot(data = RNA_data_filtLocal, aes(x = log10(p.value), y = log10(pvalue_old))) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local eQTLs Left Ventricle") +
  xlab("log10(p-values) using relatedness as covariate")+ylab("log10(pvalues) without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-30,1,5), limits = c(-30, 1))+
  scale_y_continuous(breaks = seq(-30,0,5), limits = c(-30, 0))+geom_text(aes(-28, 0 , label=paste0("R = ", round(cor.test(log10(RNA_data_filtLocal$p.value), log10(RNA_data_filtLocal$pvalue_old))$estimate, 3)), vjust=0))

e=ggplot(data = Ribo_data_filtLocal, aes(x = log10(p.value), y = log10(pvalue_old))) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local riboQTLs Left Ventricle") +
  xlab("log10(p-values) using relatedness as covariate")+ylab("log10(pvalues) without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-30,1,5), limits = c(-30, 1))+
  scale_y_continuous(breaks = seq(-30,0,5), limits = c(-30, 0))+geom_text(aes(-28, 0 , label=paste0("R = ", round(cor.test(log10(Ribo_data_filtLocal$p.value), log10(Ribo_data_filtLocal$pvalue_old))$estimate, 3)), vjust=0))

f=ggplot(data = te_data_filtLocal, aes(x = log10(p.value), y = log10(pvalue_old))) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local teQTLs Left Ventricle") +
  xlab("log10(p-values) using relatedness as covariate")+ylab("log10(pvalues) without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-20,1,2.5), limits = c(-20, 1))+
  scale_y_continuous(breaks = seq(-20,0,2.5), limits = c(-20, 0))+geom_text(aes(-19, 0 , label=paste0("R = ", round(cor.test(log10(te_data_filtLocal$p.value), log10(te_data_filtLocal$pvalue_old))$estimate, 3)), vjust=0))


grid.arrange(d,e,f,a,b,c, ncol=3, nrow=2)
dev.off()

