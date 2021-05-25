#!/usr/bin/Rscript


setwd("RiboTaper/results/pooled_all_liver/")


files<-list.files(pattern="ORFs_max",all.files=T,path="RiboTaper/results/liver/results",include.dirs=F,recursive=T,full.names=T)
files_ok<-files[-unique(c(grep(x=files,pattern="_filt"),grep(x=files,pattern="pooled_all")))]
samples<-sapply(strsplit(files_ok,"/"),"[[",9)

files_pooled<-files[grep(x=files,pattern="pooled_all")]
files_pooled<-files_pooled[-(c(grep(x=files,pattern="_filt"),grep(x=files,pattern="_hg19"),grep(x=files,pattern="pooled_all_noSHR05")))]
orfs_pooled<-read.table(files_pooled,stringsAsFactors=F,header=T,sep="\t")
orfs_pooled<-orfs_pooled[which((orfs_pooled$pct_covered_onlymulti_ribo/orfs_pooled$pct_region_covered_ribo)<0.3),]

inters_strand<-read.table("orfs_inters_str",stringsAsFactors=F,header=F,sep="\t")

list_samples<-list()
	for(i in 1:length(files_ok)){
	orfs<-read.table(files_ok[i],stringsAsFactors=F,header=T,sep="\t")
	orfs$sample<-samples[i]
	orfs[orfs$strand=="+","gen_end"]<-unlist(lapply((lapply(strsplit(orfs[orfs$strand=="+","ORF_id_gen"],"_"),"[",c(1,3))),paste,collapse="_"))
	orfs[orfs$strand=="-","gen_end"]<-unlist(lapply((lapply(strsplit(orfs[orfs$strand=="-","ORF_id_gen"],"_"),"[",c(1,2))),paste,collapse="_"))
	list_samples[[i]]<-orfs
}

all_orfs<-do.call(list_samples,what=rbind.data.frame)

orfs_pooled[orfs_pooled$strand=="+","gen_end" ]<-unlist(lapply((lapply(strsplit(orfs_pooled[orfs_pooled$strand=="+","ORF_id_gen"],"_"),"[",c(1,3))),paste,collapse="_"))
orfs_pooled[orfs_pooled$strand=="-","gen_end"]<-unlist(lapply((lapply(strsplit(orfs_pooled[orfs_pooled$strand=="-","ORF_id_gen"],"_"),"[",c(1,2))),paste,collapse="_"))
n_samples_genstop<-c()
n_samples_ORF<-c()
n_samples_gene<-c() #Val

ids<-all_orfs$sample

genpos<-all_orfs$gen_end
orfid<-all_orfs$ORF_id_gen
geneid<- all_orfs$gene_id #Val

for(j in 1:dim(orfs_pooled)[1]){
	id<-orfs_pooled[j,c("ORF_id_gen","gen_end","gene_id")] #Val, added 'gene_id'
	n_samples_genstop[j]<-length(unique(ids[genpos%in%id[2]]))
	n_samples_ORF[j]<-length(unique(ids[orfid%in%id[1]]))
	n_samples_gene[j] <- length(unique(ids[geneid%in%id[3]])) #Val
}



orfs_pooled$samples_exact_ORF<-n_samples_ORF
orfs_pooled$samples_stop_ORF<-n_samples_genstop
orfs_pooled$samples_gene<-n_samples_gene #Val
orfs_pooled$overlap_cds_samestr<-orfs_pooled$ORF_id_tr%in%inters_strand$V1
write.table(orfs_pooled,file="pooled_ORFs_table",sep="\t",row.names=F,col.names=T,quote=F)

orfs_bed<-read.table("translated_ORFs_sorted.bed",stringsAsFactors=F)
uorfs<-orfs_pooled[which(orfs_pooled$category=="uORF" & orfs_pooled$ORF_P_sites>20 & orfs_pooled$ORF_length>23),]

orfs_select<-orfs_bed[which(ids%in%uorfs$ORF_id_tr),]
write.table(uorfs,file="uORFs_selected",quote=F,col.names=T,row.names=F,sep="\t")
write.table(orfs_select,file="uorfs_selected.bed",quote=F,col.names=F,row.names=F,sep="\t")
system('less uorfs_selected.bed | sort -k1,1 -k2,2n > uORFs_selected_sorted.bed')
system('rm uorfs_selected.bed')



dorfs<-orfs_pooled[which(orfs_pooled$category=="dORF" & orfs_pooled$ORF_P_sites>20 & orfs_pooled$ORF_length>23),]

orfs_select<-orfs_bed[which(ids%in%dorfs$ORF_id_tr),]
write.table(dorfs,file="dORFs_selected",quote=F,col.names=T,row.names=F,sep="\t")
write.table(orfs_select,file="dorfs_selected.bed",quote=F,col.names=F,row.names=F,sep="\t")
system('less dorfs_selected.bed | sort -k1,1 -k2,2n > dORFs_selected_sorted.bed')
system('rm dorfs_selected.bed')


