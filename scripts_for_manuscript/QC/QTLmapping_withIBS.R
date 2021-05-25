res_dir="/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/"


load(file = paste(res_dir, "Result_tables/20180213_Lv_polyA/eqtls.RData", sep=""))
qtls_lv_rna = actual.eqtls
load(file=paste(res_dir, "Result_tables/20180213_Lv_Ribo/eqtls.RData", sep=""))
qtls_lv_ribo = actual.eqtls
load(file=paste(res_dir, "Result_tables/20180213_Lv_Ribo_RNA_resid/eqtls.RData", sep=""))
qtls_lv_te = actual.eqtls

load(file = paste(res_dir, "Result_tables/20201105_Lv_polyA_ibs/eqtls.RData", sep=""))
qtls_lv_rna_ibs = actual.eqtls
load(file=paste(res_dir, "Result_tables/20201105_Lv_Ribo_ibs/eqtls.RData", sep=""))
qtls_lv_ribo_ibs = actual.eqtls
load(file=paste(res_dir, "Result_tables/20201107_Lv_Ribo_RNA_resid_ibs/eqtls.RData", sep=""))
qtls_lv_te_ibs = actual.eqtls


qtls_lv_rna = rbind(qtls_lv_rna[["cis"]], qtls_lv_rna[["trans"]])
qtls_lv_rna_ibs = rbind(qtls_lv_rna_ibs[["cis"]], qtls_lv_rna_ibs[["trans"]])
qtls_lv_ribo = rbind(qtls_lv_ribo[["cis"]], qtls_lv_ribo[["trans"]])
qtls_lv_ribo_ibs = rbind(qtls_lv_ribo_ibs[["cis"]], qtls_lv_ribo_ibs[["trans"]])
qtls_lv_te = rbind(qtls_lv_te[["cis"]], qtls_lv_te[["trans"]])
qtls_lv_te_ibs = rbind(qtls_lv_te_ibs[["cis"]], qtls_lv_te_ibs[["trans"]])



qtls_lv_rna$key = paste(qtls_lv_rna$SNP, qtls_lv_rna$gene, sep="_")
qtls_lv_rna_ibs$key = paste(qtls_lv_rna_ibs$SNP, qtls_lv_rna_ibs$gene, sep="_")

RNA_data = data.frame(key1 = qtls_lv_rna$key, SNP = qtls_lv_rna$SNP, gene = qtls_lv_rna$gene, beta_old = qtls_lv_rna$beta, 
                      pvalue_old = qtls_lv_rna$`p-value`, beta_ibs = qtls_lv_rna_ibs[match(qtls_lv_rna$key, qtls_lv_rna_ibs$key), "beta"], 
                      pvalue_ibs = qtls_lv_rna_ibs[match(qtls_lv_rna$key, qtls_lv_rna_ibs$key), "p-value"],
                      key2 = qtls_lv_rna_ibs[match(qtls_lv_rna$key, qtls_lv_rna_ibs$key), "key"] )

lv_rna_local = read.table("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/local_eQTLs_Lv_polyA_sign.txt", sep="\t", header=T, stringsAsFactors = F)
lv_rna_distant = read.table("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/distant_eQTLs_Lv_polyA_sign.txt", sep="\t", header=T, stringsAsFactors = F)
lv_rna_local$key = paste(lv_rna_local$SNP, lv_rna_local$gene, sep="_")
lv_rna_distant$key = paste(lv_rna_distant$SNP, lv_rna_distant$gene, sep="_")

RNA_data_filtLocal = RNA_data[match(lv_rna_local$key, RNA_data$key1),]
RNA_data_filtDistant = RNA_data[match(lv_rna_distant$key, RNA_data$key1),]

cor.test(RNA_data_filtLocal$beta, RNA_data_filtLocal$beta_old) #0.9666
cor.test(RNA_data_filtLocal$p.value, RNA_data_filtLocal$pvalue_old, method = "spearman") #0.564

cor.test(RNA_data_filtDistant$beta, RNA_data_filtDistant$beta_old) # 0.9913
cor.test(RNA_data_filtDistant$p.value, RNA_data_filtDistant$pvalue_old,method = "spearman") #0.6497


qtls_lv_ribo$key = paste(qtls_lv_ribo$SNP, qtls_lv_ribo$gene, sep="_")
qtls_lv_ribo_ibs$key = paste(qtls_lv_ribo_ibs$SNP, qtls_lv_ribo_ibs$gene, sep="_")

Ribo_data = data.frame(key1 = qtls_lv_ribo$key, SNP = qtls_lv_ribo$SNP, gene_name = qtls_lv_ribo$gene_name, gene = qtls_lv_ribo$gene, beta_old = qtls_lv_ribo$beta, 
                      pvalue_old = qtls_lv_ribo$`p-value`, beta_ibs = qtls_lv_ribo_ibs[match(qtls_lv_ribo$key, qtls_lv_ribo_ibs$key), "beta"], 
                      pvalue_ibs = qtls_lv_ribo_ibs[match(qtls_lv_ribo$key, qtls_lv_ribo_ibs$key), "p-value"],
                      key2 = qtls_lv_ribo_ibs[match(qtls_lv_ribo$key, qtls_lv_ribo_ibs$key), "key"] )

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
write.table(Ribo_data_filtDistant, file = "Desktop/Ribo_data_filtDistant.txt", sep="\t", col.names=T, row.names = F, quote = F)
#save.image("Desktop/MatrixEQTL_withIBS.RData")





qtls_lv_te$key = paste(qtls_lv_te$SNP, qtls_lv_te$gene, sep="_")
qtls_lv_te_ibs$key = paste(qtls_lv_te_ibs$SNP, qtls_lv_te_ibs$gene, sep="_")

te_data = data.frame(key1 = qtls_lv_te$key, SNP = qtls_lv_te$SNP, gene_name = qtls_lv_te$gene_name, gene = qtls_lv_te$gene, beta_old = qtls_lv_te$beta, 
                       pvalue_old = qtls_lv_te$`p-value`, beta_ibs = qtls_lv_te_ibs[match(qtls_lv_te$key, qtls_lv_te_ibs$key), "beta"], 
                       pvalue_ibs = qtls_lv_te_ibs[match(qtls_lv_te$key, qtls_lv_te_ibs$key), "p-value"],
                       key2 = qtls_lv_te_ibs[match(qtls_lv_te$key, qtls_lv_te_ibs$key), "key"] )

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
write.table(te_data_filtDistant, file = "Desktop/TE_data_filtDistant.txt", sep="\t", col.names=T, row.names = F, quote = F)
save.image("Desktop/MatrixEQTL_withIBS.RData")


library(ggplot2)
pdf("Desktop/Scatter_QTLMapping_CompareIBS_noIBS.pdf")
ggplot(data = RNA_data_filtDistant, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant eQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.5))+
  scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.5))+geom_text(aes(-1.3, 1.5 , label=paste0("R = ", round(cor.test(RNA_data_filtDistant$beta, RNA_data_filtDistant$beta_old)$estimate, 3)), vjust=0))
ggplot(data = Ribo_data_filtDistant, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant riboQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.7))+
  scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.7))+geom_text(aes(-1.3, 1.7 , label=paste0("R = ", round(cor.test(Ribo_data_filtDistant$beta, Ribo_data_filtDistant$beta_old)$estimate, 3)), vjust=0))
ggplot(data = te_data_filtDistant, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Distant teQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.5))+
  scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5, 1.5))+geom_text(aes(-1.3, 1.5 , label=paste0("R = ", round(cor.test(te_data_filtDistant$beta, te_data_filtDistant$beta_old)$estimate, 3)), vjust=0))

ggplot(data = RNA_data_filtLocal, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local eQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-4,4,1), limits = c(-4,4))+
  scale_y_continuous(breaks = seq(-4,4,1), limits = c(-4,4)) +geom_text(aes(-3.5, 4 , label=paste0("R = ", round(cor.test(RNA_data_filtLocal$beta, RNA_data_filtLocal$beta_old)$estimate, 3)), vjust=0))

ggplot(data = Ribo_data_filtLocal, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local riboQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-3,5.5,1), limits = c(-3,5))+
  scale_y_continuous(breaks =  seq(-3,5,1), limits = c(-3,5)) +geom_text(aes(-2.5, 5 , label=paste0("R = ", round(cor.test(Ribo_data_filtLocal$beta, Ribo_data_filtLocal$beta_old)$estimate, 3)), vjust=0))

ggplot(data = te_data_filtLocal, aes(x = beta, y = beta_old)) + geom_point(alpha = 0.3, col = "blue") + ggtitle("Local teQTLs Left Ventricle") +
  xlab("effect size using relatedness as covariate")+ylab("effect size without covariates")+
  theme_classic()+geom_abline(intercept = 0, slope = 1)+scale_x_continuous(breaks = seq(-1.5,1.5,0.5), limits=c(-2,1.6))+
  scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits=c(-2,1.6)) +geom_text(aes(-1.7, 1.6 , label=paste0("R = ", round(cor.test(te_data_filtLocal$beta, te_data_filtLocal$beta_old)$estimate, 3)), vjust=0))

dev.off()

