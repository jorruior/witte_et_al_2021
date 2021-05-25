source("Process_CountData.R")

library(seqinr)

prefix = "path_to_data"

add_fpkm <- function(tab, fpkm_tab, colname){
  fpkm_tab_mean <- apply(fpkm_tab,1,mean) # calculate mean
  fpkm_tab_mean_df <- data.frame(gene_id=names(fpkm_tab_mean), FPKM_mean=as.numeric(fpkm_tab_mean))
  tab[,colname] <- fpkm_tab_mean_df[match(tab$gene_id, fpkm_tab_mean_df$gene_id),"FPKM_mean"] # add to table
  return(tab)
}

##main function
filterTable <- function(orfs, rna_gene_fpkm, sub, out){
  #remove unwanted columns
  rem_col=c("P_sites_sum","RNA_sites","Ribo_cov_aver","RNA_cov_aver","annotated_start","annotated_stop","ORF_reads_ribo","ORF_reads_rna","ORF_RNA_sites",
            "ORF_RNAsit_pct_in_frame","ORF_pval_multi_rna","ORF_spec_multi_ribo","ORF_spec_multi_rna","ORF_id_tr_annotated","pct_region_covered_ribo",
            "pct_covered_onlymulti_ribo","pct_region_covered_rna","pct_covered_onlymulti_rna","Method")
  orfs <- orfs[ , setdiff(colnames(orfs),rem_col)]

  #add RNA FPKM to the table
  orfs =  add_fpkm(tab=orfs, fpkm_tab=rna_gene_fpkm, colname="FPKM_mean_rna_gene")
  # filter table
  orfs = subset(orfs, samples_gene >= 10 & FPKM_mean_rna_gene>=1)
 
 orfs_filt_df <- data.frame(lapply(orfs, as.character), stringsAsFactors=FALSE)
  write.table(orfs_filt_df, file=out, sep="\t", col.names=T, row.names=F, quote=F)
  
  
  return(orfs)
}


RI_lv_all_ORFs <- read.table(file.path(prefix,"20171204_RatHeartTranslatome/results/RiboTaper/results/pooled_all_lv/pooled_ORFs_table"), sep="\t", header=T)
RI_lv_rna_gene_fpkm  <- read.table(file.path(prefix,"20171204_RatHeartTranslatome/data/count_matrices/FPKM/polyA_RI_lv_proteincodingANDlincRNA_countMatrix_FPKM.txt"), header=T)
out_path <- file.path(prefix,"20171204_RatHeartTranslatome/results/Detected_ORFs/RI_lv_all_ORFs_filtered.txt")
all_filtered_orfs = filterTable(orfs=RI_lv_all_ORFs, rna_gene_fpkm=RI_lv_rna_gene_fpkm, sub="rat", out=out_path)

lncRNA_filtered_orfs = all_filtered_orfs[which((all_filtered_orfs$category=="ncORFS"| !all_filtered_orfs$overlap_cds_samestr ) & all_filtered_orfs$annotation %in% c("stringtie","antisense","lincRNA","processed_transcript") & all_filtered_orfs$ORF_P_sites>20 & all_filtered_orfs$ORF_length>23),]

write.table(lncRNA_filtered_orfs, file = paste(prefix,"RatHeartTranslatome/results/Detected_ORFs/RI_lv_lncRNA_ORFs_filtered.txt", sep="/"), sep="\t", col.names=T, row.names=F, quote=F)



RI_li_all_ORFs <- read.table(file.path(prefix,"RatHeartTranslatome/results/RiboTaper/results/pooled_all_liver/pooled_ORFs_table"), sep="\t", header=T)
RI_li_rna_gene_fpkm  <- read.table( file.path(prefix,"RatHeartTranslatome/data/count_matrices/FPKM/polyA_RI_li_proteincodingANDlincRNA_countMatrix_FPKM.txt"), header=T)
out_path <- file.path(prefix,"RatHeartTranslatome/results/Detected_ORFs/RI_li_all_ORFs_filtered.txt")
all_filtered_orfs_li = filterTable(orfs=RI_li_all_ORFs, rna_gene_fpkm=RI_li_rna_gene_fpkm, sub="rat", out=out_path)

lncRNA_filtered_orfs_li = all_filtered_orfs_li[which((all_filtered_orfs_li$category=="ncORFS"| !all_filtered_orfs_li$overlap_cds_samestr ) & all_filtered_orfs_li$annotation %in% c("stringtie","antisense","lincRNA","processed_transcript") & all_filtered_orfs_li$ORF_P_sites>20 & all_filtered_orfs_li$ORF_length>23),]
write.table(lncRNA_filtered_orfs_li, file = paste(prefix,"RatHeartTranslatome/results/Detected_ORFs/RI_li_lncRNA_ORFs_filtered.txt", sep="/"), sep="\t", col.names=T, row.names=F, quote=F)
