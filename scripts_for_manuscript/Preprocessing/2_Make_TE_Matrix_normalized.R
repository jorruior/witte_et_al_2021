library(DESeq2)


data_dir = "RatHeartTranslatome/data/count_matrices/"

#Load count matrices (raw data)

#RI lv
raw_counts_ribo_RI_lv = read.table(file=paste(data_dir, "ribo_RI_lv_proteincodingANDlincRNA_countMatrix_raw.txt", sep=""), header=T, sep="\t")
raw_counts_polyA_RI_lv = read.table(file=paste(data_dir, "polyA_RI_lv_proteincodingANDlincRNA_countMatrix_raw.txt", sep=""), header=T, sep="\t")

#RI li
raw_counts_ribo_RI_li = read.table(file=paste(data_dir, "ribo_RI_li_proteincodingANDlincRNA_countMatrix_raw.txt", sep=""), header=T, sep="\t")
raw_counts_polyA_RI_li = read.table(file=paste(data_dir, "polyA_RI_li_proteincodingANDlincRNA_countMatrix_raw.txt", sep=""), header=T, sep="\t")

#congenics lv

raw_counts_ribo_cong_lv = read.table(file=paste(data_dir, "ribo_cong_lv_proteincodingANDlincRNA_countMatrix_raw.txt", sep=""), header=T, sep="\t")
raw_counts_polyA_cong_lv = read.table(file=paste(data_dir, "polyA_cong_lv_proteincodingANDlincRNA_countMatrix_raw.txt", sep=""), header=T, sep="\t")

#function to make normalized translational efficiency matrix
TE_normalization = function(rna, ribo, pattern)
{
	pool_sizeFactor <- estimateSizeFactorsForMatrix(cbind(rna,ribo))
	
	mrna_sizeFactor <- pool_sizeFactor[1:ncol(rna)]
	rpf_sizeFactor <- pool_sizeFactor[(ncol(rna)+1):(ncol(rna)+ncol(ribo))]
	polyA_norm = t(apply(rna, 1, function(x) {x/mrna_sizeFactor}))
	ribo_norm = t(apply(ribo, 1, function(x) {x/rpf_sizeFactor}))
	TE = ribo_norm/polyA_norm

	write.table(polyA_norm, file=paste("RatHeartTranslatome/data/count_matrices/deseq_norm/", pattern, "_RNA_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
	write.table(ribo_norm, file=paste("RatHeartTranslatome/data/count_matrices/deseq_norm/", pattern, "_Ribo_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
	write.table(TE, file=paste("RatHeartTranslatome/data/count_matrices/TE/", pattern, "_TE_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
	
	return(TE)

}


TE_polyA_RI_lv = TE_normalization(raw_counts_polyA_RI_lv, raw_counts_ribo_RI_lv, "polyA_RI_lv")
TE_polyA_RI_li = TE_normalization(raw_counts_polyA_RI_li, raw_counts_ribo_RI_li, "polyA_RI_li")
TE_polyA_cong_lv = TE_normalization(raw_counts_polyA_cong_lv, raw_counts_ribo_cong_lv, "polyA_cong_lv")










