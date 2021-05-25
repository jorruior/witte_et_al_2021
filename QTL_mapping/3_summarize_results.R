input_dir = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/final_table/"
#lv
load(paste(input_dir, "Lv_polyA_tables.RData", sep=""))
lv_polyA_local = eQTL_cis.sign
lv_polyA_distant = eQTL_trans.sign

load(paste(input_dir, "Lv_Ribo_tables.RData", sep=""))
lv_ribo_local = eQTL_cis.sign
lv_ribo_distant = eQTL_trans.sign

load(paste(input_dir, "Lv_RiboRNAresid_tables.RData", sep=""))
lv_te_local = eQTL_cis.sign
lv_te_distant = eQTL_trans.sign

all_sets_lv = list(lv_polyA_local = unique(lv_polyA_local$gene), lv_polyA_distant = unique(lv_polyA_distant$gene), lv_ribo_local = unique(lv_ribo_local$gene), lv_ribo_distant = unique(lv_ribo_distant$gene), lv_te_local = unique(lv_te_local$gene), lv_te_distant = unique(lv_te_distant$gene))
all_sets_lv_noTE = list(lv_polyA_local = unique(lv_polyA_local$gene), lv_polyA_distant = unique(lv_polyA_distant$gene), lv_ribo_local = unique(lv_ribo_local$gene), lv_ribo_distant = unique(lv_ribo_distant$gene))


#liver
load(paste(input_dir, "Liver_polyA_tables.RData", sep=""))
li_polyA_local = eQTL_cis.sign
li_polyA_distant = eQTL_trans.sign

load(paste(input_dir, "Liver_Ribo_tables.RData", sep=""))
li_ribo_local = eQTL_cis.sign
li_ribo_distant = eQTL_trans.sign

load(paste(input_dir, "Liver_RiboRNAresid_tables.RData", sep=""))
li_te_local = eQTL_cis.sign
li_te_distant = eQTL_trans.sign

all_sets_li = list(li_polyA_local = unique(li_polyA_local$gene), li_polyA_distant = unique(li_polyA_distant$gene), li_ribo_local = unique(li_ribo_local$gene), li_ribo_distant = unique(li_ribo_distant$gene), li_te_local = unique(li_te_local$gene), li_te_distant = unique(li_te_distant$gene))

all_sets_li_noTE = list(li_polyA_local = unique(li_polyA_local$gene), li_polyA_distant = unique(li_polyA_distant$gene), li_ribo_local = unique(li_ribo_local$gene), li_ribo_distant = unique(li_ribo_distant$gene))


#same results
library(gplots)
ItemsList_lv <- venn(all_sets_lv_noTE, show.plot=FALSE)
ItemsList_li <- venn(all_sets_li_noTE, show.plot=FALSE)

save.image(paste(input_dir, "All_tables.RData", sep=""))
