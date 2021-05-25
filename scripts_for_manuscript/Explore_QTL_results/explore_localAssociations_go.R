source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/SeqQTL_oldCluster.R")

library(GenomicRanges)
load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/permutations/All_tables_permutations.RData")

load(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/filtered_genesets_ribotaper_fpkm.RData")

#go enrichment of locally regulated gene sets

library(gProfileR)
set_base_url("http://biit.cs.ut.ee/gprofiler_archive2/r1477_e82_eg29/web")

local_eQTL_heart_go = gprofiler(query =  unique(lv_polyA_local$gene), custom_bg = genesset_lv, organism = "rnorvegicus")
local_riboQTL_heart_go = gprofiler(query =  unique(lv_ribo_local$gene), custom_bg = genesset_lv, organism = "rnorvegicus")
local_eQTL_liver_go = gprofiler(query =  unique(li_polyA_local$gene), custom_bg = genesset_li, organism = "rnorvegicus")
local_riboQTL_liver_go = gprofiler(query =  unique(li_ribo_local$gene), custom_bg = genesset_li, organism = "rnorvegicus")

## not very conclusive.


#compare results with parental findings
load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/DESeq_parentals/deseq2_results_parentals_li.RData")
load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/DESeq_parentals/deseq2_results_parentals_lv.RData")

sum(unique(lv_polyA_local$gene) %in% rownames(diff_res_polyA_lv)) #702 
sum(unique(lv_ribo_local$gene) %in% rownames(diff_res_ribo_lv)) #309
sum(unique(li_polyA_local$gene) %in% rownames(diff_res_polyA_li)) #365 
sum(unique(li_ribo_local$gene) %in% rownames(diff_res_ribo_li)) #383

#all genes with local QTLs found in parentals as differential as well
#same effect trend?

sum(unique(lv_polyA_distant$gene) %in% rownames(diff_res_polyA_lv)) #32 
sum(unique(lv_ribo_distant$gene) %in% rownames(diff_res_ribo_lv)) #45
sum(unique(li_polyA_distant$gene) %in% rownames(diff_res_polyA_li)) #12 
sum(unique(li_ribo_distant$gene) %in% rownames(diff_res_ribo_li)) #0
#all genes with distant QTLs found in parentals as differential as well

#use table from sebas and ele

differential_genes_ele <- read.table("/home/fkreuch/work/Ribsome_Profiling/data/slope_all_new.txt", header=T, sep="\t")

sum(unique(lv_polyA_local$gene) %in% as.character(differential_genes_ele[which((differential_genes_ele$Groups %in% c("Both", "RNA exclusive")) & (differential_genes_ele$Tissue == "LV sep")), "X"])) #170/702
sum(unique(li_polyA_local$gene) %in% as.character(differential_genes_ele[which((differential_genes_ele$Groups %in% c("Both", "RNA exclusive")) & (differential_genes_ele$Tissue == "Liver Seq")), "X"])) #125/365
sum(unique(lv_ribo_local$gene) %in% as.character(differential_genes_ele[which((differential_genes_ele$Groups %in% c("Both", "RPF exclusive")) & (differential_genes_ele$Tissue == "LV sep")), "X"])) #136/309
sum(unique(li_ribo_local$gene) %in% as.character(differential_genes_ele[which((differential_genes_ele$Groups %in% c("Both", "RPF exclusive")) & (differential_genes_ele$Tissue == "Liver Seq")), "X"])) #142/383


sum(unique(lv_polyA_distant$gene) %in% as.character(differential_genes_ele[which((differential_genes_ele$Groups %in% c("Both", "RNA exclusive")) & (differential_genes_ele$Tissue == "LV sep")), "X"])) #2/32
sum(unique(li_polyA_distant$gene) %in% as.character(differential_genes_ele[which((differential_genes_ele$Groups %in% c("Both", "RNA exclusive")) & (differential_genes_ele$Tissue == "Liver Seq")), "X"])) #4/12
sum(unique(lv_ribo_distant$gene) %in% as.character(differential_genes_ele[which((differential_genes_ele$Groups %in% c("Both", "RPF exclusive")) & (differential_genes_ele$Tissue == "LV sep")), "X"])) #5/45
sum(unique(li_ribo_distant$gene) %in% as.character(differential_genes_ele[which((differential_genes_ele$Groups %in% c("Both", "RPF exclusive")) & (differential_genes_ele$Tissue == "Liver Seq")), "X"])) #0/0




distant_eQTL_heart_go = gprofiler(query =  unique(lv_polyA_distant$gene), custom_bg = genesset_lv, organism = "rnorvegicus")
distant_riboQTL_heart_go = gprofiler(query =  unique(lv_ribo_distant$gene), custom_bg = genesset_lv, organism = "rnorvegicus")
distant_eQTL_liver_go = gprofiler(query =  unique(li_polyA_distant$gene), custom_bg = genesset_li, organism = "rnorvegicus")
distant_riboQTL_liver_go = gprofiler(query =  unique(li_ribo_distant$gene), custom_bg = genesset_li, organism = "rnorvegicus")
#nothings


