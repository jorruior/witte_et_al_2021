
library(eQTLpipeline)
library(data.table)
library(GenomicRanges)
library(grid)

source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/eqtl_trans_kinship.R")
source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/SeqQTL_oldCluster.R")
load("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment/Input_MatrixEQTL/Liver_uORF/Liver_Ribo_uORF_input.RData")
uORF_ribo_liver = pheno

load("/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment/Input_MatrixEQTL/Lv_uORF/Lv_Ribo_uORF_input.RData")
uORF_ribo_lv = pheno



load(file="/Volumes/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/preprocessed_data.RData")
sample_names = c("BXH02", "BXH03", "BXH05", "HXB01", "HXB02", "HXB03", "HXB04", "HXB05", "HXB07", "HXB10", "HXB13", "HXB15", "HXB17", "HXB18", "HXB20", "HXB21", "HXB22", "HXB23", "HXB24", "HXB25", "HXB27", "HXB29", "HXB31", "BXH06", "BXH08", "BXH09", "BXH10", "BXH11", "BXH12", "BXH13") 

#Phenotype Input file for eQTL Mapping
log2_ribo_lv = log2_ribo_lv[,match(sample_names, colnames(log2_ribo_lv))]
colnames(log2_ribo_li) = gsub("_4", "", colnames(log2_ribo_li))
log2_ribo_li = log2_ribo_li[,match(sample_names, colnames(log2_ribo_li))]

log2_te_lv = resid(lm(as.matrix(log2_ribo_lv)~as.matrix(log2_polyA_lv)))
log2_te_li = resid(lm(as.matrix(log2_ribo_li)~as.matrix(log2_polyA_li)))



load("/Volumes//huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment/20180227_image.RData")
load("/Volumes//huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/permutations/All_tables_permutations.RData")

polyA_lv_with_uORF = data.frame(lv_polyA_local[which(lv_polyA_local$gene %in% all_uORFs_lv$gene_id),])
ribo_lv_with_uORF = data.frame(lv_ribo_local[which(lv_ribo_local$gene %in% all_uORFs_lv$gene_id),])
te_lv_with_uORF = data.frame(lv_te_local[which(lv_te_local$gene %in% all_uORFs_lv$gene_id),])


polyA_li_with_uORF = data.frame(li_polyA_local[which(li_polyA_local$gene %in% all_uORFs_li$gene_id),])
ribo_li_with_uORF = data.frame(li_ribo_local[which(li_ribo_local$gene %in% all_uORFs_li$gene_id),])
te_li_with_uORF = data.frame(li_te_local[which(li_te_local$gene %in% all_uORFs_li$gene_id),])


wilcoxon_uORF = function(uORF_expression, QTL_results, pattern, genotypes_alleles,all_uORFs,mORF_expression)
{
	diffExprResult = NULL	
	for(i in 1:nrow(QTL_results))
	{
		SDP = QTL_results[i,]$SNP		
		gene_id = QTL_results[i,]$gene
		genotypes = genotypes_alleles[SDP,]

		uORFs = all_uORFs[which(all_uORFs$gene_id == gene_id),"ORF_id_tr"]

		expression = uORF_expression[na.omit(match(uORFs, rownames(uORF_expression))),]
		
		
		for (uORF in rownames(expression))
		{
			expression2 = expression[uORF,which(genotypes != "1")]			
			status = factor(genotypes[which(genotypes != "1")], levels = c("0","2"))
			FC = mean(as.numeric(expression[uORF,which(status == "2")])) / mean(as.numeric(expression[uORF,which(status == "0")]))
			mORF_FC = mean(as.numeric(mORF_expression[gene_id,which(status == "2")])) / mean(as.numeric(mORF_expression[gene_id,which(status == "0")]))		
	
			p = wilcox.test(as.numeric(expression2[uORF,]) ~ status)$p.value
			
			diffExprResult = rbind(diffExprResult, data.frame(gene_id = gene_id, SNP = SDP,uORF_id = uORF, pvalue = p, uORF_FC = FC, mORF_FC = mORF_FC))
			
		}
		
		

	}
	diffExprResult$FDR = p.adjust(diffExprResult$pvalue, method ="BH")
	
	return(diffExprResult)
}

uORF_link_mORFQTL_ribo_lv = wilcoxon_uORF(uORF_ribo_lv, ribo_lv_with_uORF, "riboQTL_lv", genotypes_alleles,all_uORFs_lv, log2_ribo_lv)
uORF_link_mORFQTL_ribo_lv.sign = uORF_link_mORFQTL_ribo_lv[which(uORF_link_mORFQTL_ribo_lv$FDR < 0.05),] #38
uORF_link_mORFQTL_TE_lv = wilcoxon_uORF(uORF_ribo_lv, te_lv_with_uORF, "teQTL_lv", genotypes_alleles,all_uORFs_lv, log2_te_lv)
uORF_link_mORFQTL_TE_lv.sign = uORF_link_mORFQTL_TE_lv[which(uORF_link_mORFQTL_TE_lv$FDR < 0.05),] #10
uORF_link_mORFQTL_ribo_li = wilcoxon_uORF(uORF_ribo_liver, ribo_li_with_uORF, "riboQTL_li", genotypes_alleles, all_uORFs_li,log2_ribo_li)
uORF_link_mORFQTL_ribo_li.sign = uORF_link_mORFQTL_ribo_li[which(uORF_link_mORFQTL_ribo_li$FDR < 0.05),]#38
uORF_link_mORFQTL_TE_li = wilcoxon_uORF(uORF_ribo_liver, te_li_with_uORF, "teQTL_li", genotypes_alleles,all_uORFs_li, log2_te_li) #nothing significant


uORF_link_mORFQTL_ribo_lv.sign[which(uORF_link_mORFQTL_ribo_lv.sign$mORF_FC < 1 & uORF_link_mORFQTL_ribo_lv.sign$uORF_FC > 1),]
uORF_link_mORFQTL_ribo_lv.sign[which(uORF_link_mORFQTL_ribo_lv.sign$uORF_FC < 1 & uORF_link_mORFQTL_ribo_lv.sign$mORF_FC > 1),]

uORF_link_mORFQTL_TE_lv.sign[which(uORF_link_mORFQTL_TE_lv.sign$mORF_FC < 1 & uORF_link_mORFQTL_TE_lv.sign$uORF_FC > 1),]
uORF_link_mORFQTL_TE_lv.sign[which(uORF_link_mORFQTL_TE_lv.sign$uORF_FC < 1 & uORF_link_mORFQTL_TE_lv.sign$mORF_FC > 1),]

uORF_link_mORFQTL_ribo_li.sign[which(uORF_link_mORFQTL_ribo_li.sign$mORF_FC < 1 & uORF_link_mORFQTL_ribo_li.sign$uORF_FC > 1),]
uORF_link_mORFQTL_ribo_li.sign[which(uORF_link_mORFQTL_ribo_li.sign$uORF_FC < 1 & uORF_link_mORFQTL_ribo_li.sign$mORF_FC > 1),]

library(cowplot)
library(gridExtra)

pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment/uORF_with_TEQTL_lv_opposite.pdf", width=10)

for(i in 1:nrow(uORF_link_mORFQTL_TE_lv.sign))
{

	hit = uORF_link_mORFQTL_TE_lv.sign[i,]
	SDP =hit$SNP
	gene_name = hit$gene_name
	genotype = genotypes_alleles[SDP,]
	gene = hit$gene
	gene1 = data.frame(expression = c(t(log2_ribo_lv[gene,]), log2_te_lv[gene,], t(log2_polyA_lv[gene,]), t(uORF_ribo_lv[gene,])), dataset = c(rep("ribo", 30), rep("te", 30), rep("polyA", 30), rep("uORF", 30)), genotype = rep(genotype, 4))
	
	gene1 = gene1[which(gene1$genotype %in% c("0", "2")),]
	colors = c("black", "red")

	plot1 = qplot(as.factor(genotype), expression, geom="boxplot", data=gene1[which(gene1$dataset %in% c("polyA", "ribo")),], main = paste("Expression of ", hit$ORF_id,  " ~ Strain", sep=""), col=genotype, xlab="", ylab="normalized expression values") + facet_grid(.~dataset) + theme_bw() + scale_colour_manual(values=colors)

	plot2 = qplot(as.factor(genotype), expression, geom="boxplot", data=gene1[which(gene1$dataset %in% c("te")),], col=genotype, xlab="", ylab="TE") + facet_grid(.~dataset) + theme_bw() + scale_colour_manual(values=colors)
	
	plot3 = qplot(as.factor(genotype), expression, geom="boxplot", data=gene1[which(gene1$dataset %in% c("uORF")),], col=genotype, xlab="", ylab="uORF ribo") + facet_grid(.~dataset) + theme_bw() + scale_colour_manual(values=colors)

	
t = plot_grid(plot1, plot2, plot3, align = "h", nrow = 1, rel_widths = c(0.4,0.2,0.2))
 print(t)

}

dev.off()

pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment/uORF_with_riboQTL_lv_opposite.pdf", width=10)

for(i in 1:nrow(uORF_link_mORFQTL_ribo_lv.sign))
{

	hit = uORF_link_mORFQTL_ribo_lv.sign[i,]
	SDP =hit$SNP
	gene_name = hit$gene_name
	genotype = genotypes_alleles[SDP,]
	gene = hit$gene
	gene1 = data.frame(expression = c(t(log2_ribo_lv[gene,]), log2_te_lv[gene,], t(log2_polyA_lv[gene,]), t(uORF_ribo_lv[gene,])), dataset = c(rep("ribo", 30), rep("te", 30), rep("polyA", 30), rep("uORF", 30)), genotype = rep(genotype, 4))
	
	gene1 = gene1[which(gene1$genotype %in% c("0", "2")),]
	colors = c("black", "red")

	plot1 = qplot(as.factor(genotype), expression, geom="boxplot", data=gene1[which(gene1$dataset %in% c("polyA", "ribo")),], main = paste("Expression of ", hit$ORF_id,  " ~ Strain", sep=""), col=genotype, xlab="", ylab="normalized expression values") + facet_grid(.~dataset) + theme_bw() + scale_colour_manual(values=colors)

	plot2 = qplot(as.factor(genotype), expression, geom="boxplot", data=gene1[which(gene1$dataset %in% c("te")),], col=genotype, xlab="", ylab="TE") + facet_grid(.~dataset) + theme_bw() + scale_colour_manual(values=colors)
	
	plot3 = qplot(as.factor(genotype), expression, geom="boxplot", data=gene1[which(gene1$dataset %in% c("uORF")),], col=genotype, xlab="", ylab="uORF ribo") + facet_grid(.~dataset) + theme_bw() + scale_colour_manual(values=colors)

	
t = plot_grid(plot1, plot2, plot3, align = "h", nrow = 1, rel_widths = c(0.4,0.2,0.2))
 print(t)

}

dev.off()



pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment/uORF_with_riboQTL_li_opposite.pdf", width=10)

for(i in 1:nrow(uORF_link_mORFQTL_ribo_li.sign))
{

	hit = uORF_link_mORFQTL_ribo_li.sign[i,]
	SDP =hit$SNP
	gene_name = hit$gene_name
	genotype = genotypes_alleles[SDP,]
	gene = hit$gene
	gene1 = data.frame(expression = c(t(log2_ribo_li[gene,]), log2_te_li[gene,], t(log2_polyA_li[gene,]), t(uORF_ribo_liver[gene,])), dataset = c(rep("ribo", 30), rep("te", 30), rep("polyA", 30), rep("uORF", 30)), genotype = rep(genotype, 4))
	
	gene1 = gene1[which(gene1$genotype %in% c("0", "2")),]
	colors = c("black", "red")

	plot1 = qplot(as.factor(genotype), expression, geom="boxplot", data=gene1[which(gene1$dataset %in% c("polyA", "ribo")),], main = paste("Expression of ", hit$ORF_id,  " ~ Strain", sep=""), col=genotype, xlab="", ylab="normalized expression values") + facet_grid(.~dataset) + theme_bw() + scale_colour_manual(values=colors)

	plot2 = qplot(as.factor(genotype), expression, geom="boxplot", data=gene1[which(gene1$dataset %in% c("te")),], col=genotype, xlab="", ylab="TE") + facet_grid(.~dataset) + theme_bw() + scale_colour_manual(values=colors)
	
	plot3 = qplot(as.factor(genotype), expression, geom="boxplot", data=gene1[which(gene1$dataset %in% c("uORF")),], col=genotype, xlab="", ylab="uORF ribo") + facet_grid(.~dataset) + theme_bw() + scale_colour_manual(values=colors)

	
t = plot_grid(plot1, plot2, plot3, align = "h", nrow = 1, rel_widths = c(0.4,0.2,0.2))
 print(t)

}

dev.off()


