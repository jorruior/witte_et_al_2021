##FC Ribo vs. Length
library(RSvgDevice)

load(file = "/VOLUMES/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/DESeq_congenics/result_deseq.RData")
load(file="/VOLUMES/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/filtered_genesets_ribotaper_fpkm.RData")

diff_res_ribo_lv$FC = 2^(diff_res_ribo_lv$log2FoldChange)
diff_res_polyA_lv$FC = 2^(diff_res_polyA_lv$log2FoldChange)

Ribotaper="/VOLUMES/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/"
ORFs_lv = read.table(paste(Ribotaper, "RI_lv_all_ORFs_filtered.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)
filtered_ORFs_lv = ORFs_lv[which(ORFs_lv$category %in% c("ORFs_ccds", "nonccds_coding_ORFs", "ncORFS") ),]

FC_CDS_length = data.frame(gene_id = genesset_lv, Ribo_FC = diff_res_ribo_lv[genesset_lv,]$FC, RNA_FC = diff_res_polyA_lv[genesset_lv,]$FC, CDS_length = sapply(genesset_lv, function(x) max(filtered_ORFs_lv[which(filtered_ORFs_lv$gene_id == x),"ORF_length"])))


ribo_deseq_cong_lv = read.table(file = "/Volumes//huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/deseq_norm/polyA_cong_lv_Ribo_deseqNorm_RNAANDRIBO.txt", sep="\t", header=T)
polya_deseq_cong_lv = read.table(file = "/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/deseq_norm/polyA_cong_lv_RNA_deseqNorm_RNAANDRIBO.txt", sep="\t", header=T)

residuals=resid(lm(as.matrix(log2(ribo_deseq_cong_lv[genesset_lv,]+1))~as.matrix(log2(polya_deseq_cong_lv[genesset_lv,]+1))))
residuals = residuals[match(rownames(FC_CDS_length), rownames(residuals)),]
residuals_mean = data.frame(BN_like = apply(residuals[,1:5], 1, mean), SHR_like = apply(residuals[,6:10], 1, mean))
FC_residual = residuals_mean$SHR_like / residuals_mean$BN_like

FC_CDS_length$FC_residual = FC_residual

#FC_CDS_length$FC_residual = abs(log2(FC_CDS_length$Ribo_FC)) - abs(log2(FC_CDS_length$RNA_FC))

library(reshape2)
library(ggplot2)

FC_CDS_length.melted = melt(FC_CDS_length, id = c("gene_id", "CDS_length"))


library(smatr)
library(scales)

#pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/DESeq_congenics/length_effect_onlyCongenics.pdf", width=9, height=5)
devSVG("/VOLUMES/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/DESeq_congenics/length_effect_onlyCongenics_new.svg",width=15, height=5)

par(mfrow=c(1,3))
test1 <- sma(FC_CDS_length$CDS_length~FC_CDS_length$Ribo_FC, log="xy")

plot(10^(as.numeric(test1$data[,"FC_CDS_length$Ribo_FC"])), 10^(as.numeric(test1$data[,"FC_CDS_length$CDS_length"])), ylab="CDS_length", xlab="riboFC(SHR/BN)", main="Congenics riboFC",  type="p", log="xy", xlim=c(0.2,5), pch=19)
plot(test1, pch=20, col="red", type="l", add= TRUE)
abline(0,0,v=1, col="blue")
abline(0,0,v=1.2, col="green", lty="dashed")
abline(0,0,v=1/1.2, col="green", lty="dashed")

test2 <- sma(FC_CDS_length$CDS_length~FC_CDS_length$RNA_FC, log="xy")

plot(10^(as.numeric(test2$data[,"FC_CDS_length$RNA_FC"])), 10^(as.numeric(test2$data[,"FC_CDS_length$CDS_length"])), ylab="CDS_length", xlab="polyAFC(SHR/BN)", main="Congenics polyAFC",  type="p", log="xy", xlim=c(0.2,5), pch=19)
abline(0,0,v=1, col="blue")
abline(0,0,v=1.2, col="green", lty="dashed")
abline(0,0,v=1/1.2, col="green", lty="dashed")

test3 <- sma(FC_CDS_length$CDS_length~FC_CDS_length$FC_residual, log="xy")

plot(10^(as.numeric(test3$data[,"FC_CDS_length$FC_residual"])), 10^(as.numeric(test3$data[,"FC_CDS_length$CDS_length"])), ylab="CDS_length", xlab="residual FC(SHR/BN)", main="Congenics residual FC",  type="p", log="xy", xlim=c(0.2,5), pch=19)
plot(test3, pch=20, col="red", type="l", add= TRUE)

abline(0,0,v=1, col="blue")
abline(0,0,v=1.2, col="green", lty="dashed")
abline(0,0,v=1/1.2, col="green", lty="dashed")

dev.off()


