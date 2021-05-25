input_dir = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/permutations/"
#lv

lv_polyA_local = read.table(paste(input_dir, "Lv_polyA_esnps_cis.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)
lv_polyA_distant = read.table(paste(input_dir, "Lv_polyA_esnps_trans.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)

lv_polyA_local = lv_polyA_local[which(lv_polyA_local$FDR_new <= 0.1 & lv_polyA_local$empirical.p <= 0.0015),]
lv_polyA_distant = lv_polyA_distant[which(lv_polyA_distant$FDR_new <= 0.2 & lv_polyA_distant$empirical.p <= 0.0015),]

lv_ribo_local = read.table(paste(input_dir, "Lv_ribo_esnps_cis.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)
lv_ribo_distant = read.table(paste(input_dir, "Lv_ribo_esnps_trans.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)


lv_ribo_local = lv_ribo_local[which(lv_ribo_local$FDR_new <= 0.1 & lv_ribo_local$empirical.p <= 0.0015),]
lv_ribo_distant = lv_ribo_distant[which(lv_ribo_distant$FDR_new <= 0.2 & lv_ribo_distant$empirical.p <= 0.0015),]


lv_te_local = read.table(paste(input_dir, "Lv_RiboRNAresid_esnps_cis.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)
lv_te_distant = read.table(paste(input_dir, "Lv_RiboRNAresid_esnps_trans.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)

lv_te_local = lv_te_local[which(lv_te_local$FDR_new <= 0.1 & lv_te_local$empirical.p <= 0.0015),]
lv_te_distant = lv_te_distant[which(lv_te_distant$FDR_new <= 0.2 & lv_te_distant$empirical.p <= 0.0015),]


all_sets_lv_local = list(lv_polyA_local = unique(lv_polyA_local$gene),  lv_ribo_local = unique(lv_ribo_local$gene), lv_te_local = unique(lv_te_local$gene))
all_sets_lv_distant = list(lv_polyA_distant = unique(lv_polyA_distant$gene),  lv_ribo_distant = unique(lv_ribo_distant$gene), lv_te_distant = unique(lv_te_distant$gene))


all_sets_lv_noTE = list(lv_polyA_local = unique(lv_polyA_local$gene), lv_polyA_distant = unique(lv_polyA_distant$gene), lv_ribo_local = unique(lv_ribo_local$gene), lv_ribo_distant = unique(lv_ribo_distant$gene))


#liver

li_polyA_local = read.table(paste(input_dir, "Liver_polyA_esnps_cis.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)
li_polyA_distant = read.table(paste(input_dir, "Liver_polyA_esnps_trans.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)

li_polyA_local = li_polyA_local[which(li_polyA_local$FDR_new <= 0.1 & li_polyA_local$empirical.p <= 0.0015),]
li_polyA_distant = li_polyA_distant[which(li_polyA_distant$FDR_new <= 0.2 & li_polyA_distant$empirical.p <= 0.0015),]

li_ribo_local = read.table(paste(input_dir, "Liver_Ribo_esnps_cis.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)
li_ribo_distant = read.table(paste(input_dir, "Liver_Ribo_esnps_trans.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)


li_ribo_local = li_ribo_local[which(li_ribo_local$FDR_new <= 0.1 & li_ribo_local$empirical.p <= 0.0015),]
li_ribo_distant = li_ribo_distant[which(li_ribo_distant$FDR_new <= 0.2 & li_ribo_distant$empirical.p <= 0.0015),]


li_te_local = read.table(paste(input_dir, "Liver_RiboRNAresid_esnps_cis.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)
li_te_distant = read.table(paste(input_dir, "Liver_RiboRNAresid_esnps_trans.txt", sep=""), sep="\t", header=T, stringsAsFactor=F)


li_te_local = li_te_local[which(li_te_local$FDR_new <= 0.1 & li_te_local$empirical.p <= 0.0015),]
li_te_distant = li_te_distant[which(li_te_distant$FDR_new <= 0.2 & li_te_distant$empirical.p <= 0.0015),]




all_sets_li_local = list(li_polyA_local = unique(li_polyA_local$gene),  li_ribo_local = unique(li_ribo_local$gene), li_te_local = unique(li_te_local$gene)) 
all_sets_li_distant = list(li_polyA_distant = unique(li_polyA_distant$gene),  li_ribo_distant = unique(li_ribo_distant$gene), li_te_distant = unique(li_te_distant$gene))

all_sets_li_noTE = list(li_polyA_local = unique(li_polyA_local$gene), li_polyA_distant = unique(li_polyA_distant$gene), li_ribo_local = unique(li_ribo_local$gene), li_ribo_distant = unique(li_ribo_distant$gene))



#same results
library(gplots)
ItemsList_lv <- venn(all_sets_lv_noTE, show.plot=FALSE)
ItemsList_li <- venn(all_sets_li_noTE, show.plot=FALSE)

ItemsList_lv_te_local<- venn(all_sets_lv_local, show.plot=TRUE)
ItemsList_lv_te_distant<- venn(all_sets_lv_distant, show.plot=TRUE)
ItemsList_li_te_local<- venn(all_sets_li_local, show.plot=TRUE)
ItemsList_li_te_distant<- venn(all_sets_li_distant, show.plot=TRUE)

save.image(paste(input_dir, "All_tables_permutations.RData", sep=""))


write.table(lv_polyA_local, file = paste(input_dir, "Lv_local_eQTL.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
write.table(lv_ribo_local, file = paste(input_dir, "Lv_local_riboQTL.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
write.table(lv_te_local, file = paste(input_dir, "Lv_local_teQTL.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
write.table(lv_polyA_distant, file = paste(input_dir, "Lv_distant_eQTL.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
write.table(lv_ribo_distant, file = paste(input_dir, "Lv_distant_riboQTL.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
write.table(lv_te_distant, file = paste(input_dir, "Lv_distant_teQTL.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)

write.table(li_polyA_local, file = paste(input_dir, "Liver_local_eQTL.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
write.table(li_ribo_local, file = paste(input_dir, "Liver_local_riboQTL.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
write.table(li_te_local, file = paste(input_dir, "Liver_local_teQTL.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
write.table(li_polyA_distant, file = paste(input_dir, "Liver_distant_eQTL.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
write.table(li_ribo_distant, file = paste(input_dir, "Liver_distant_riboQTL.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
write.table(li_te_distant, file = paste(input_dir, "Liver_distant_teQTL.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=F)
