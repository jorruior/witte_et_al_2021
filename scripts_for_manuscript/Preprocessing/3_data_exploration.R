library(eQTLpipeline)
library(preprocessCore)
library(GenomicRanges)
library(DESeq2)
library(reshape)
library(ggplot2)
library(plyr)
library(scales)

source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/SeqQTL_oldCluster.R")
source("~/work/HS_heart/R/functions/vif_function.R")


#get ensembl gene annotation
res_dir="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/"
count_matrices="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/"
Ribotaper="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/"

genes_names = gtf2GR("/data/huebner2/DATASETS/GENOMES/Ensembl_Genomes/Ens82/rn6/annotation/Rattus_norvegicus.Rnor_6.0.82.gtf", c("gene_id", "gene_name", "gene_biotype"))
genes_names = genes_names[values(genes_names)[,"type"] == "gene"]
genes_names = genes_names[which(seqnames(genes_names) %in% c(1:20, "X")),]



polyA_fpkm_lv = read.table(paste(count_matrices, "FPKM/polyA_RI_lv_proteincodingANDlincRNA_countMatrix_FPKM.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)
polyA_deseq_lv = read.table(paste(count_matrices, "deseq_norm/polyA_RI_lv_RNA_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)

polyA_fpkm_li = read.table(paste(count_matrices, "FPKM/polyA_RI_li_proteincodingANDlincRNA_countMatrix_FPKM.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)
polyA_deseq_li = read.table(paste(count_matrices, "deseq_norm/polyA_RI_li_RNA_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)


ribo_fpkm_lv = read.table(paste(count_matrices, "FPKM/ribo_RI_lv_proteincodingANDlincRNA_countMatrix_FPKM.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)
ribo_deseq_lv = read.table(paste(count_matrices, "deseq_norm/polyA_RI_lv_Ribo_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)

ribo_fpkm_li = read.table(paste(count_matrices, "FPKM/ribo_RI_li_proteincodingANDlincRNA_countMatrix_FPKM.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)
ribo_deseq_li = read.table(paste(count_matrices, "deseq_norm/polyA_RI_li_Ribo_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)



a_polyA_lv = apply(polyA_fpkm_lv,1,function(x){expressed=mean(x)>1})
a_polyA_liver = apply(polyA_fpkm_li,1,function(x){expressed=mean(x)>1})


filtered_fpkm_polyA_lv= polyA_fpkm_lv[a_polyA_lv,]
filtered_fpkm_polyA_li= polyA_fpkm_li[a_polyA_liver,]

pdf(file=paste(res_dir, "expressed_genes_FPKM1.pdf", sep=""))

hist(as.matrix(log10(polyA_fpkm_lv)),breaks=50,col=rgb(0,0,1,0.5),main="",xlab="Log10 FPKM Genes RI Strains LV PolyA" )
hist(as.matrix(log10(filtered_fpkm_polyA_lv)),breaks=50,add=T,col=rgb(0,0,0,0.5))

hist(as.matrix(log10(polyA_fpkm_li)),breaks=100,col=rgb(0,0,1,0.5),main="",xlab="Log10 FPKM Genes RI Strains Liver PolyA" )
hist(as.matrix(log10(filtered_fpkm_polyA_li)),breaks=100,add=T,col=rgb(0,0,0,0.5))

dev.off()

ORFs_lv = read.table(paste(Ribotaper, "RI_lv_all_ORFs_filtered.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)
ORFs_li = read.table(paste(Ribotaper, "RI_li_all_ORFs_filtered.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)

filtered_ORFs_lv = ORFs_lv[which(ORFs_lv$category %in% c("ORFs_ccds", "nonccds_coding_ORFs", "ncORFS") ),]
filtered_ORFs_li = ORFs_li[which(ORFs_li$category %in% c("ORFs_ccds", "nonccds_coding_ORFs", "ncORFS") ),]
genesset_lv = unique(filtered_ORFs_lv$gene_id)
genesset_li = unique(filtered_ORFs_li$gene_id)



filtered_deseq_polyA_lv= polyA_deseq_lv[genesset_lv,]
filtered_deseq_ribo_lv= ribo_deseq_lv[genesset_lv,]
filtered_deseq_polyA_li= polyA_deseq_li[genesset_li,]
filtered_deseq_ribo_li= ribo_deseq_li[genesset_li,]

save(filtered_deseq_polyA_lv, filtered_deseq_polyA_li, filtered_deseq_ribo_lv, filtered_deseq_ribo_li, file = paste(res_dir, "DESeq_norm_filtered_matrices_FPKM1.RData", sep=""))

make_plots_biotypes_filtered_genessets = function(genesset, out_name)
{
	table_biotypes <- as.data.frame(table(genes_names[which(genes_names$gene_id %in% genesset),]$gene_biotype))
	table_biotypes$group = "miscRNA"
	table_biotypes[which(table_biotypes$Var1 == "protein_coding"),]$group = "protein_coding"
	table_biotypes[which(table_biotypes$Var1 %in% c("antisense", "lincRNA", "processed_transcript")),]$group = "lncRNA"
	table_biotypes[which(table_biotypes$Var1 %in% c("processed_pseudogene", "pseudogene", "transcribed_processed_pseudogene", "unprocessed_pseudogene")),]$group = "pseudogene"

	forPlot = data.frame(row.names = table_biotypes$Var1, numberGenes = table_biotypes$Freq, group = table_biotypes$group)
	forPlot.melted = melt(forPlot)
	forPlot.melted$percentage = 0
	forPlot.melted$bigGroup = "plot2"
	forPlot.melted[which(forPlot.melted$group %in% c("protein_coding")),]$bigGroup = "plot1"
	forPlot.melted = rbind(forPlot.melted, forPlot.melted[which(forPlot.melted$bigGroup == "plot2"),])
	forPlot.melted$bigGroup = as.character(forPlot.melted$bigGroup)
	forPlot.melted$bigGroup[1:9] = "plot1"
	forPlot.melted[which(forPlot.melted$bigGroup == "plot1"),]$percentage = 100*(forPlot.melted[which(forPlot.melted$bigGroup == "plot1"),]$value / sum(forPlot.melted[which(forPlot.melted$bigGroup == "plot1"),]$value))
	forPlot.melted[which(forPlot.melted$bigGroup == "plot2"),]$percentage = 100*(forPlot.melted[which(forPlot.melted$bigGroup == "plot2"),]$value / sum(forPlot.melted[which(forPlot.melted$bigGroup == "plot2"),]$value))
	forPlot.melted$group = factor(forPlot.melted$group, levels = c("protein_coding", "lncRNA", "pseudogene", "miscRNA"))

	pdf(paste(res_dir, out_name, sep=""))
	p = ggplot() +
	  geom_bar(aes(y = percentage, x = bigGroup, fill = group), data = forPlot.melted, stat="identity") +
	  theme(legend.position="bottom", legend.direction="horizontal",
		legend.title = element_blank()) +
	  scale_y_continuous(labels = dollar_format(suffix = "%", prefix = "")) +
	  labs(x="", y="Percentage") +
	  ggtitle("Number of translated genes (%)") +
	   theme(axis.line = element_line(size=1, colour = "black"),
		panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
		panel.border = element_blank(), panel.background = element_blank())
		print(p)
	dev.off()
}

make_plots_biotypes_filtered_genessets(genesset_lv, "filtered_biotypes_lv.pdf")
make_plots_biotypes_filtered_genessets(genesset_li, "filtered_biotypes_li.pdf")


correlation_test = function(counts, prefix)
{
	test = cor(counts, method="pearson")^2
	library(ggplot2)
	library(reshape)
	library(gridExtra)
	melted_test = melt(test)	
	diag(test) = NA
	means_samples_corr= apply(test, 1, function(x){mean(x,na.rm=TRUE)})
	print(means_samples_corr[means_samples_corr < 0.9])
	pdf(file=paste(res_dir, "sample_correlation_", prefix, ".pdf", sep=""), width=24, height=8)
	plot1 = ggplot(melted_test, aes(X1, X2)) +    geom_tile(aes(fill = value)) + 
    geom_text(aes(fill = melted_test$value, label = round(melted_test$value, 5))) +
    scale_fill_gradient(low = "white", high = "red", limits=c(0.92,1))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
	print(plot1)
	#pairs(log10(counts))
	
	dev.off()
	return(test)
}

corr_polyA_lv = correlation_test(filtered_deseq_polyA_lv, "polyA_lv")
corr_polyA_li = correlation_test(filtered_deseq_polyA_li, "polyA_li")
corr_ribo_lv = correlation_test(filtered_deseq_ribo_lv, "ribo_lv")
corr_ribo_li = correlation_test(filtered_deseq_ribo_li, "ribo_li")

### lv ribo HXB03 and HXB31 have r² < 0.9 but above 0.85


write.table(corr_polyA_lv, file = paste(res_dir, "sample_correlation_polyA_lv_pearsonr2.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
write.table(corr_polyA_li, file = paste(res_dir, "sample_correlation_polyA_li_pearsonr2.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
write.table(corr_ribo_lv, file = paste(res_dir, "sample_correlation_ribo_lv_pearsonr2.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
write.table(corr_ribo_li, file = paste(res_dir, "sample_correlation_ribo_li_pearsonr2.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)



###replicates

ribo_raw_li = read.table(paste(count_matrices, "raw/ribo_RI_li_proteincodingANDlincRNA_countMatrix_raw.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)

ribo_norm_li_replicates = ribo_raw_li[,match(c("BXH12_2_li_CDS", "BXH12_3_li_CDS", "BXH12_4_li_CDS", "BXH13_2_li_CDS", "BXH13_3_li_CDS", "BXH13_4_li_CDS"), colnames(ribo_raw_li))]
covar_ribo_norm_li_replicates = data.frame(tissue = rep("liver", 6), row.names = colnames(ribo_norm_li_replicates))
summary_ribo_li_rep = DESeqDataSetFromMatrix(countData=ribo_norm_li_replicates, colData=covar_ribo_norm_li_replicates,  design=~1)
dds_ribo_li_rep <- estimateSizeFactors(summary_ribo_li_rep)
normCounts_ribo_li_rep <- counts(dds_ribo_li_rep, normalized=TRUE)

normCounts_ribo_li_rep = normCounts_ribo_li_rep[genesset_li,]


save.image(file = paste(res_dir, "replicates_data.RData", sep=""))


png(file=paste(res_dir, "sample_correlation_ribo_deseqNorm_replicates_liver.png", sep=""), width = 6.25,
  height    = 6.25,
  units     = "in",
  res       = 1200,
  pointsize = 4)

par(mar=c(3,3,3,3),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.main  = 2
)
	counts=normCounts_ribo_li_rep
	test = cor(counts, method="pearson")
	# make layout for plot layout
	laymat<-diag(1:ncol(counts)) #histograms
	laymat[upper.tri(laymat)]<-(ncol(counts)+1):(ncol(counts)+((ncol(counts)^2 - ncol(counts))/2)) #correlations
	laymat[lower.tri(laymat)]<-((ncol(counts)+((ncol(counts)^2 - ncol(counts))/2))+1):(ncol(counts)+((ncol(counts)^2 - ncol(counts))/2)+((ncol(counts)^2 - ncol(counts))/2)) #heatmaps
	layout(laymat) #define layout using laymat
# Draw histograms
for(i in 1:ncol(counts)) 
{  #hist(log10(counts[,i]+1),main=colnames(counts)[i], breaks=50, , xlab="", ylab="")
plot.new()	
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="white")
text(x=0.5,y=0.5,labels=colnames(counts)[i],cex=2)
}
##corr
hmcol = colorRampPalette(c("darksalmon", "darkred"))(n = 15)
cor_levels = levels(cut(test, 10))
codes = matrix(match(cut(test, 10), cor_levels), nrow=ncol(counts))

for(i in 2:ncol(counts))
  for(j in 1:(i-1)){
  	plot.new()
	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = hmcol[codes[i,j]])
	text(x=0.5,y=0.5,labels=round(test[i,j], 5),cex=3)
    } 
# Smoothscatter
for(i in 1:(ncol(counts)-1))
  for(j in (i+1):ncol(counts)){
    plot(log10(counts[,i]+1), log10(counts[,j]+1), xlab="", ylab="")
   
    }
dev.off()





###cross experimental correlation for each sample
library(reshape)

genesset_both = intersect(genesset_lv, genesset_li)

filtered2_deseq_polyA_lv= as.data.frame(polyA_deseq_lv[genesset_both,])
filtered2_deseq_ribo_lv= ribo_deseq_lv[genesset_both,]
filtered2_deseq_polyA_li= polyA_deseq_li[genesset_both,]
filtered2_deseq_ribo_li= ribo_deseq_li[genesset_both,]


all_datasets = data.frame(genesset_both, melt(filtered2_deseq_polyA_lv), melt(filtered2_deseq_ribo_lv)[,2], melt(filtered2_deseq_polyA_li)[,2], melt(filtered2_deseq_ribo_li)[,2])

colnames(all_datasets) = c("genes", "samples", "polyA_lv", "ribo_lv", "polyA_li", "ribo_li")
all_datasets$samples = gsub("_TV__3_CDS", "", all_datasets$samples)
all_datasets$samples = gsub("_TV__4_CDS", "", all_datasets$samples)

pdf(file=paste(res_dir, "cross_experiment_correlation.pdf", sep=""))
library(gplots)
hmcol = colorRampPalette(c("darkred", "white", "darkblue"))(n = 10)
hclustfunc <- function(x) hclust(x, method="complete")
for(sample in unique(all_datasets$samples))
{
heatmap.2(cor(as.matrix(all_datasets[which(all_datasets$samples == sample), -c(1,2)]), method="pearson")^2,col=hmcol, trace="none", symbreaks=T, margin=c(8,8), dendrogram="none", main=sample, hclust=hclustfunc, cellnote=as.matrix(round(cor(as.matrix(all_datasets[which(all_datasets$samples == sample), -c(1,2)])), 3)))
}
dev.off()

#####bis hier durch



covariates_liver_pola = read.table(file="~/work/new_data_RI_translational_efficiency/Franzi/data/current/covariates/covariates_liver_pola_min.txt", sep="\t", header=T, stringsAsFactors=F)
covariates_lv_pola = read.table(file="~/work/new_data_RI_translational_efficiency/Franzi/data/current/covariates/covariates_lv_pola_min.txt", sep="\t", header=T, stringsAsFactors=F)
covariates_lv_ribo = read.table("~/work/RI_translational_efficiency/data/20151113/covariates_RPF/RPF_covariates_lv.txt", sep="\t", header=T, stringsAsFactors = F)
covariates_liver_ribo = read.table("~/work/RI_translational_efficiency/data/20151113/covariates_RPF/RPF_covariates_liver.txt", sep="\t", header=T, stringsAsFactors = F)



PCA = function(log2_counts, covariates)
{
  library(reshape)
  pca <- prcomp(t(log2_counts))
  predict_pca <- predict(pca)

  pvalue_matrix = apply(predict_pca, 2, function(pc) {sapply(covariates[1:ncol(covariates)], function(covar) {
    avail = !is.na(covar)
    m0 = lm(pc[avail] ~ 1)
    m1 = lm(pc[avail] ~ covar[avail])
    return(anova(m1, m0)[2,"Pr(>F)"])
  })
  })
  pvalue_matrix[is.na(pvalue_matrix)] = NA
  melted_pvalue_matrix = melt(pvalue_matrix)
  colnames(melted_pvalue_matrix) = c("covar", "PC", "pvalue")
  melted_pvalue_matrix.plot = na.omit(melted_pvalue_matrix)
  
  return(melted_pvalue_matrix.plot)
  
}




#polya

covariates_lv_pola = covariates_lv_pola[match(colnames(filtered_deseq_polyA_lv), gsub("pol", "CDS", covariates_lv_pola$sample_name)),]
covariates_lv_pola.min = data.frame(row.names = covariates_lv_pola$sample_name, model.matrix(covariates_lv_pola$sample_name ~ covariates_lv_pola$RNA_isolation_date), model.matrix(covariates_lv_pola$sample_name ~ covariates_lv_pola$X1st_seq_date), covariates_lv_pola$RIN_score)[,-1]
pca_covar_lv_polya = PCA(filtered_deseq_polyA_lv, covariates_lv_pola.min)
summary(prcomp(t(filtered_deseq_polyA_lv))) #1st PC explains 46.7% of variance
pca_covar_lv_polya[pca_covar_lv_polya$pvalue < 0.01,] #no covariates 


covariates_li_pola = covariates_liver_pola[match(colnames(filtered_deseq_polyA_li), gsub("pol", "CDS", covariates_liver_pola$sample_name)),]
covariates_li_pola.min = data.frame(row.names = covariates_li_pola$sample_name, model.matrix(covariates_li_pola$sample_name ~ covariates_li_pola$RNA_isolation_date), model.matrix(covariates_li_pola$sample_name ~ covariates_li_pola$X1st_seq_date), covariates_li_pola$RIN_score)[,-1]
pca_covar_li_polya = PCA(filtered_deseq_polyA_li, covariates_li_pola.min)
summary(prcomp(t(filtered_deseq_polyA_li))) #1st PC explains 80.9% of variance
pca_covar_li_polya[pca_covar_li_polya$pvalue < 0.01,] # no covariates

#ribo

covariates_lv_ribo = covariates_lv_ribo[match(gsub("_lv_CDS", "", colnames(filtered_deseq_ribo_lv)), covariates_lv_ribo$sample_ID),]
covariates_lv_ribo.min = data.frame(row.names = covariates_lv_ribo$sample_ID, model.matrix(covariates_lv_ribo$sample_ID ~ covariates_lv_ribo$produced_by), model.matrix(covariates_lv_ribo$sample_ID ~ covariates_lv_ribo$date_of_powdering), model.matrix(covariates_lv_ribo$sample_ID ~ covariates_lv_ribo$date_of_PCR), model.matrix(covariates_lv_ribo$sample_ID ~ covariates_lv_ribo$no_of_PCR_cycles),  covariates_lv_ribo$library_.concentration_ng_uL)[,-1]
pca_covar_lv_ribo = PCA(filtered_deseq_ribo_lv, covariates_lv_ribo.min) 
summary(prcomp(t(filtered_deseq_ribo_lv))) #1st PC explains 61.6% of variance
pca_covar_lv_ribo[pca_covar_lv_ribo$pvalue < 0.01,] # no covariates


covariates_li_ribo = covariates_liver_ribo[match(gsub("_li_CDS", "", colnames(filtered_deseq_ribo_li)), covariates_liver_ribo$sample_ID),]
covariates_li_ribo.min = data.frame(row.names = covariates_li_ribo$sample_ID, model.matrix(covariates_li_ribo$sample_ID ~ covariates_li_ribo$date_of_powdering), model.matrix(covariates_li_ribo$sample_ID ~ covariates_li_ribo$date_of_PCR), model.matrix(covariates_li_ribo$sample_ID ~ covariates_li_ribo$no_of_PCR_cycles),  covariates_li_ribo$library_.concentration_ng_uL)[,-1]
pca_covar_li_ribo = PCA(filtered_deseq_ribo_li, covariates_li_ribo.min) 
summary(prcomp(t(filtered_deseq_ribo_li))) #1st PC explains 81.5% of variance
pca_covar_li_ribo[pca_covar_li_ribo$pvalue < 0.01,] #no covariates



log_transform=function(qnorm_counts)
{
	a = qnorm_counts
	a[a=="0"]=1
	qnorm_counts[qnorm_counts=="0"]=min(a)
	log2_qnorm=log2(qnorm_counts)
	rownames(log2_qnorm) = rownames(qnorm_counts)
	colnames(log2_qnorm) = colnames(qnorm_counts)

	return(log2_qnorm)
}

log2_ribo_li = log_transform(filtered_deseq_ribo_li)
log2_ribo_lv = log_transform(filtered_deseq_ribo_lv)
log2_polyA_li = log_transform(filtered_deseq_polyA_li)
log2_polyA_lv = log_transform(filtered_deseq_polyA_lv)



colnames(log2_ribo_li) = colnames(log2_ribo_lv) = colnames(log2_polyA_li) = colnames(log2_polyA_lv) = gsub("_lv_CDS", "", colnames(log2_ribo_lv))


pdf(file=paste(res_dir, "normalized_counts_FPKM1.pdf", sep=""), width=10)
par(mar=c(10.1,4.1,4.1,4.1))
boxplot(log2_ribo_li, main = "Normalized Ribo Liver counts", las=2)
boxplot(log2_ribo_lv, main = "Normalized Ribo Heart counts", las=2)
boxplot(log2_polyA_li, main = "Normalized PolyA Liver counts", las=2)
boxplot(log2_polyA_lv, main = "Normalized PolyA Heart counts", las=2)
dev.off()


sample_names = c("BXH02", "BXH03", "BXH05", "HXB01", "HXB02", "HXB03", "HXB04", "HXB05", "HXB07", "HXB10", "HXB13", "HXB15", "HXB17", "HXB18", "HXB20", "HXB21", "HXB22", "HXB23", "HXB24", "HXB25", "HXB27", "HXB29", "HXB31", "BXH06", "BXH08", "BXH09", "BXH10", "BXH11", "BXH12", "BXH13") 

log2_ribo_li = log2_ribo_li[,match(sample_names, colnames(log2_ribo_li))]
log2_ribo_lv = log2_ribo_lv[,match(sample_names, colnames(log2_ribo_lv))]
log2_polyA_li = log2_polyA_li[,match(sample_names, colnames(log2_polyA_li))]
log2_polyA_lv = log2_polyA_lv[,match(sample_names, colnames(log2_polyA_lv))]



save(log2_ribo_li, log2_ribo_lv, log2_polyA_li, log2_polyA_lv, genes_names, file=paste(res_dir, "preprocessed_data.RData", sep=""))

#kick out snoRNA and pseudogenes for FDR/QTL mapping

filtered_ORFs_lv = ORFs_lv[which(ORFs_lv$category %in% c("ORFs_ccds", "nonccds_coding_ORFs", "ncORFS") & ORFs_lv$annotation %in% c("antisense", "protein_coding", "lincRNA", "processed_transcript")),]
filtered_ORFs_li = ORFs_li[which(ORFs_li$category %in% c("ORFs_ccds", "nonccds_coding_ORFs", "ncORFS") & ORFs_li$annotation %in% c("antisense", "protein_coding", "lincRNA", "processed_transcript")),]

genesset_lv = unique(filtered_ORFs_lv$gene_id)
genesset_li = unique(filtered_ORFs_li$gene_id)

save(genesset_lv, genesset_li, file=paste(res_dir, "filtered_genesets_ribotaper_fpkm.RData", sep=""))


