res_dir="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/"

load(file=paste(res_dir, "preprocessed_data.RData", sep=""))


log2_polyA_lv.mean = apply(log2_polyA_lv, 1, function(x) mean(x))
log2_ribo_lv.mean = apply(log2_ribo_lv, 1, function(x) mean(x))

r_square_lv = cor(log2_polyA_lv.mean, log2_ribo_lv.mean, method = "pearson")^2

plot(log2_polyA_lv.mean, log2_ribo_lv.mean)

write.table(data.frame(polyA = log2_polyA_lv.mean, ribo=log2_ribo_lv.mean), file = paste(res_dir, "mean_expression_ribo_polyA_lv.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)

log2_polyA_li.mean = apply(log2_polyA_li, 1, function(x) mean(x))
log2_ribo_li.mean = apply(log2_ribo_li, 1, function(x) mean(x))

r_square_li = cor(log2_polyA_li.mean, log2_ribo_li.mean, method = "pearson")^2

write.table(data.frame(polyA = log2_polyA_li.mean, ribo=log2_ribo_li.mean), file = paste(res_dir, "mean_expression_ribo_polyA_li.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)




###compare total vs ribo vs. polyA
total_lv = read.table(file=paste("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/deseq_norm/", "total_RI_lv", "_RNA_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)
ribo_lv = read.table(file=paste("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/deseq_norm/", "total_RI_lv", "_Ribo_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)

log2_total_lv.mean = apply(log2(total_lv+1), 1, function(x) mean(x))
log2_riboT_lv.mean = apply(log2(ribo_lv+1), 1, function(x) mean(x))
r_total_lv = cor(log2_total_lv.mean, log2_riboT_lv.mean, method = "pearson")


polyA_lv = read.table(file=paste("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/deseq_norm/", "polyA_RI_lv", "_RNA_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)
ribo_lv = read.table(file=paste("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/deseq_norm/", "polyA_RI_lv", "_Ribo_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)

log2_polyA_lv.mean = apply(log2(polyA_lv+1), 1, function(x) mean(x))
log2_riboP_lv.mean = apply(log2(ribo_lv+1), 1, function(x) mean(x))
r_polyA_lv = cor(log2_polyA_lv.mean, log2_riboP_lv.mean, method = "pearson")

plot(log2_polyA_lv.mean, log2_riboP_lv.mean)


write.table(data.frame(polyA = log2_polyA_lv.mean, riboP=log2_riboP_lv.mean, total = log2_total_lv.mean, riboT=log2_riboP_lv.mean), file = paste(res_dir, "mean_expression_ribo_polyA_total_lv.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)



###compare total vs ribo vs. polyA
total_li = read.table(file=paste("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/deseq_norm/", "total_RI_li", "_RNA_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)
ribo_li = read.table(file=paste("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/deseq_norm/", "total_RI_li", "_Ribo_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)

log2_total_li.mean = apply(log2(total_li+1), 1, function(x) mean(x))
log2_riboT_li.mean = apply(log2(ribo_li+1), 1, function(x) mean(x))
r_total_li = cor(log2_total_li.mean, log2_riboT_li.mean, method = "pearson")


polyA_li = read.table(file=paste("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/deseq_norm/", "polyA_RI_li", "_RNA_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)
ribo_li = read.table(file=paste("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/deseq_norm/", "polyA_RI_li", "_Ribo_deseqNorm_RNAANDRIBO.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)

log2_polyA_li.mean = apply(log2(polyA_li+1), 1, function(x) mean(x))
log2_riboP_li.mean = apply(log2(ribo_li+1), 1, function(x) mean(x))
r_polyA_li = cor(log2_polyA_li.mean, log2_riboP_li.mean, method = "pearson")
plot(log2_polyA_li.mean, log2_riboP_li.mean)

write.table(data.frame(polyA = log2_polyA_li.mean, riboP=log2_riboP_li.mean, total = log2_total_li.mean, riboT=log2_riboP_li.mean), file = paste(res_dir, "mean_expression_ribo_polyA_total_li.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)


