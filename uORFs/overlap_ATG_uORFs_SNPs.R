library(GenomicRanges)

load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_assessment/20180227_image.RData")

local_sdps = read.table("/home/fkreuch/work/RI_translational_efficiency/data/20160318/RN6_SDPs/HXBBXH-SDP_RN60_LDB-local.txt", sep="\t", comment.char = "#", header=F, stringsAsFactors=F, colClasses="character", check.names=F)
colnames(local_sdps) = c("local_ID", "chr", "start", "end", "Marker_RN60_CSV", "Marker_count")

snps = read.table("/home/fkreuch/work/RI_translational_efficiency/data/20160318/RN6_SDPs/", sep="\t", comment.char = "#", header=F, stringsAsFactors=F, colClasses="character", check.names=F)

snps.GR= GRanges(seqnames=snps$V2, ranges = IRanges(start=as.numeric(as.character(snps$V3)), end =as.numeric(as.character(snps$V3))+1), id = snps$V1, ref=snps$V4, snp=snps$V5)

library(stringr)
all_uORFs_lv.filt$chr = str_split_fixed(all_uORFs_lv.filt$ORF_id_gen, "_", 3)[,1]
all_uORFs_lv.filt$first_pos = as.numeric(str_split_fixed(all_uORFs_lv.filt$ORF_id_gen, "_", 3)[,2])
all_uORFs_lv.filt$last_pos = as.numeric(str_split_fixed(all_uORFs_lv.filt$ORF_id_gen, "_", 3)[,3])
all_uORFs_lv.filt[which(all_uORFs_lv.filt$chr == "X"),"chr"] = 21

atg_positions_lv.GR = GRanges(seqnames = all_uORFs_lv.filt$chr, ranges = IRanges(start = all_uORFs_lv.filt$first_pos, end = all_uORFs_lv.filt$last_pos), strand = all_uORFs_lv.filt$strand, ORF_id_tr = all_uORFs_lv.filt$ORF_id_tr)

test_start = start(atg_positions_lv.GR)
test_end = end(atg_positions_lv.GR)

start(atg_positions_lv.GR[which(strand(atg_positions_lv.GR)=="-"),]) = test_end[which(strand(atg_positions_lv.GR)=="-")]-2
end(atg_positions_lv.GR[which(strand(atg_positions_lv.GR)=="-"),]) = test_end[which(strand(atg_positions_lv.GR)=="-")]+2

start(atg_positions_lv.GR[which(strand(atg_positions_lv.GR)=="+"),]) = test_start[which(strand(atg_positions_lv.GR)=="+")]-2
end(atg_positions_lv.GR[which(strand(atg_positions_lv.GR)=="+"),]) = test_start[which(strand(atg_positions_lv.GR)=="+")]+2


subsetByOverlaps(atg_positions_lv.GR, snps.GR)
subsetByOverlaps(snps.GR,atg_positions_lv.GR)
#two uORFs with SNPs in ATG 
#ENSRNOT00000026297_526_631 - Nudt21
#ENSRNOT00000084045_91_121 - Zdhhc4



all_uORFs_li.filt$chr = str_split_fixed(all_uORFs_li.filt$ORF_id_gen, "_", 3)[,1]
all_uORFs_li.filt$first_pos = as.numeric(str_split_fixed(all_uORFs_li.filt$ORF_id_gen, "_", 3)[,2])
all_uORFs_li.filt$last_pos = as.numeric(str_split_fixed(all_uORFs_li.filt$ORF_id_gen, "_", 3)[,3])
all_uORFs_li.filt[which(all_uORFs_li.filt$chr == "X"),"chr"] = 21

atg_positions_li.GR = GRanges(seqnames = all_uORFs_li.filt$chr, ranges = IRanges(start = all_uORFs_li.filt$first_pos, end = all_uORFs_li.filt$last_pos), strand = all_uORFs_li.filt$strand, ORF_id_tr = all_uORFs_li.filt$ORF_id_tr)

test_start = start(atg_positions_li.GR)
test_end = end(atg_positions_li.GR)

start(atg_positions_li.GR[which(strand(atg_positions_li.GR)=="-"),]) = test_end[which(strand(atg_positions_li.GR)=="-")]-2
end(atg_positions_li.GR[which(strand(atg_positions_li.GR)=="-"),]) = test_end[which(strand(atg_positions_li.GR)=="-")]+2

start(atg_positions_li.GR[which(strand(atg_positions_li.GR)=="+"),]) = test_start[which(strand(atg_positions_li.GR)=="+")]-2
end(atg_positions_li.GR[which(strand(atg_positions_li.GR)=="+"),]) = test_start[which(strand(atg_positions_li.GR)=="+")]+2


subsetByOverlaps(atg_positions_li.GR, snps.GR)




> genotypes_alleles["SDPG_19010912369",]
BXH02 BXH03 BXH05 HXB01 HXB02 HXB03 HXB04 HXB05 HXB07 HXB10 HXB13 HXB15 HXB17 
  "0"   "0"   "0"   "2"   "0"   "2"   "2"   "0"   "0"   "0"   "2"   "0"   "2" 
HXB18 HXB20 HXB21 HXB22 HXB23 HXB24 HXB25 HXB27 HXB29 HXB31 BXH06 BXH08 BXH09 
  "0"   "2"   "2"   "2"   "0"   "0"   "2"   "2"   "0"   "0"   "2"   "2"   "2" 
BXH10 BXH11 BXH12 BXH13 
  "0"   "2"   "0"   "0" 
> genotypes_alleles["SDPG_12013285235",]
BXH02 BXH03 BXH05 HXB01 HXB02 HXB03 HXB04 HXB05 HXB07 HXB10 HXB13 HXB15 HXB17 
  "0"   "0"   "2"   "0"   "2"   "1"   "2"   "2"   "0"   "0"   "2"   "0"   "1" 
HXB18 HXB20 HXB21 HXB22 HXB23 HXB24 HXB25 HXB27 HXB29 HXB31 BXH06 BXH08 BXH09 
  "0"   "0"   "2"   "0"   "2"   "2"   "2"   "0"   "2"   "0"   "2"   "0"   "0" 
BXH10 BXH11 BXH12 BXH13 
  "2"   "2"   "2"   "2" 

