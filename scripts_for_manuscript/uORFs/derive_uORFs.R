library(eQTLpipeline)
library(data.table)
library(GenomicRanges)
library(grid)
library(preprocessCore)
library(MASS)
source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/SeqQTL_oldCluster.R")
library(biomaRt)


library(stringr)
all_ORFs_lv = read.delim("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/RI_lv_all_ORFs_filtered.txt", sep="\t", header=T, stringsAsFactors=F)
all_ORFs_li = read.delim("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/RI_li_all_ORFs_filtered.txt", sep="\t", header=T, stringsAsFactors=F)

load(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/data_exploration/filtered_genesets_ribotaper_fpkm.RData")
all_ORFs_lv = all_ORFs_lv[which(all_ORFs_lv$gene_id %in% genesset_lv),]
all_ORFs_li = all_ORFs_li[which(all_ORFs_li$gene_id %in% genesset_li),]

indpendent_uORF_lv = all_ORFs_lv[which(all_ORFs_lv$category == "uORF" & all_ORFs_lv$samples_stop_ORF >= 10),] #n=795
indpendent_uORF_li = all_ORFs_li[which(all_ORFs_li$category == "uORF" & all_ORFs_li$samples_stop_ORF >= 10),] #n=929

genes_names = gtf2GR("/data/huebner2/DATASETS/GENOMES/Ensembl_Genomes/Ens82/rn6/annotation/Rattus_norvegicus.Rnor_6.0.82.gtf", c("gene_id", "gene_name", "gene_biotype", "transcript_id"))
test = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "rnorvegicus_gene_ensembl")

derive_overlapping_uORFs = function(all_ORFs, tissue, genes_names, test)
{
	all_ORFs$gen_start = NA
	all_ORFs[which(all_ORFs$strand == "+"),]$gen_start = paste(str_split_fixed(all_ORFs[which(all_ORFs$strand == "+"),]$ORF_id_gen, "_", 3)[,1], str_split_fixed(all_ORFs[which(all_ORFs$strand == "+"),]$ORF_id_gen, "_", 3)[,2], sep="_")
	all_ORFs[which(all_ORFs$strand == "-"),]$gen_start = paste(str_split_fixed(all_ORFs[which(all_ORFs$strand == "-"),]$ORF_id_gen, "_", 3)[,1], str_split_fixed(all_ORFs[which(all_ORFs$strand == "-"),]$ORF_id_gen, "_", 3)[,3], sep="_")
	unique_stops = names(table(all_ORFs$gen_end)[table(all_ORFs$gen_end)==1])
	unique_start = names(table(all_ORFs$gen_start)[table(all_ORFs$gen_start)==1])

	all_transcript = unique(genes_names$transcript_id)
	
	appris_tag = getBM(attributes=c('ensembl_transcript_id', 'transcript_appris'),  filters = 'ensembl_transcript_id',  values = all_transcript, mart = test)
	appris_tag.filtered = unique(appris_tag[which(appris_tag$transcript_appris != ""), "ensembl_transcript_id"])
	genes_names_utr = genes_names[which(genes_names$type == "five_prime_utr"), ]
	genes_names_cds = genes_names[which(genes_names$type == "CDS"), ]
	genes_names_stop = genes_names[which(genes_names$type == "stop_codon"), ]
	start(genes_names_stop) = start(genes_names_stop) - 2
	end(genes_names_stop) = end(genes_names_stop) + 2
	genes_names_start = genes_names[which(genes_names$type == "start_codon"), ]

	all_ORFs.GR = GRanges(seqnames = str_split_fixed(all_ORFs$ORF_id_gen, "_", 3)[,1], ranges = IRanges(start = as.numeric(str_split_fixed(all_ORFs$ORF_id_gen, "_", 3)[,2]), end = as.numeric(str_split_fixed(all_ORFs$ORF_id_gen, "_", 3)[,3])), strand = all_ORFs$strand, ORF_id_tr = all_ORFs$ORF_id_tr)

	overlaps_stop_ORFs = subsetByOverlaps(all_ORFs.GR, genes_names_stop)
	ORFs_stop_filter = all_ORFs[-which(all_ORFs$ORF_id_tr %in% overlaps_stop_ORFs$ORF_id_tr),]

	stopfiltered_ORFs.GR = all_ORFs.GR[-which(all_ORFs.GR$ORF_id_tr %in% overlaps_stop_ORFs$ORF_id_tr),]

	overlaps_start_stop_ORFs = subsetByOverlaps(stopfiltered_ORFs.GR, genes_names_start, type="start")
	ORFs_start_stop_filter = ORFs_stop_filter[-which(ORFs_stop_filter$ORF_id_tr %in% overlaps_start_stop_ORFs$ORF_id_tr),]
	stopfiltered_ORFs2.GR = stopfiltered_ORFs.GR[-which(stopfiltered_ORFs.GR$ORF_id_tr %in% overlaps_start_stop_ORFs$ORF_id_tr),]
	overlaps_5utr_ORFs = subsetByOverlaps(stopfiltered_ORFs2.GR, genes_names_utr)
	stopfiltered_5UTR_ORFs.GR = stopfiltered_ORFs2.GR[which(stopfiltered_ORFs2.GR$ORF_id_tr %in% overlaps_5utr_ORFs$ORF_id_tr),]
	utr5_transcripts = unique(subsetByOverlaps(genes_names_utr,stopfiltered_ORFs2.GR)[,"transcript_id"])
	utr5_transcripts$appris = utr5_transcripts$transcript_id %in% appris_tag.filtered
	utr5_transcripts2 = utr5_transcripts[which(utr5_transcripts$appris == "TRUE"),]
	genes_names_cds.filt = genes_names_cds[which(genes_names_cds$transcript_id %in% utr5_transcripts$transcript_id),]
	overlaps_cds_ORFs = subsetByOverlaps(stopfiltered_5UTR_ORFs.GR, genes_names_cds.filt)

	potential_uORFs = all_ORFs[which(all_ORFs$ORF_id_tr %in% overlaps_cds_ORFs$ORF_id_tr),]
	potential_uORFs = potential_uORFs[which(potential_uORFs$category != "uORF"),]
	potential_uORFs = potential_uORFs[which(potential_uORFs$gen_end %in% unique_stops),]
	potential_uORFs = potential_uORFs[which(potential_uORFs$gen_start %in% unique_start),]

	write.table(potential_uORFs, file= paste("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/RI_", tissue, "_overlappinguORFs.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F) 


	return(potential_uORFs)

}


overlapping_uORFs_lv = derive_overlapping_uORFs(all_ORFs_lv, "lv", genes_names, test) #n=330 (326 genes)
overlapping_uORFs_li = derive_overlapping_uORFs(all_ORFs_li, "li", genes_names, test) #n=330 (327 genes)
#222 overlapping uORFs shared between tissues


all_uORFs_lv = data.frame(rbind(indpendent_uORF_lv, overlapping_uORFs_lv[,1:25]), uORF_type = c(rep("independent", nrow(indpendent_uORF_lv)), rep("overlapping", nrow(overlapping_uORFs_lv))))
all_uORFs_li = data.frame(rbind(indpendent_uORF_li, overlapping_uORFs_li[,1:25]), uORF_type = c(rep("independent", nrow(indpendent_uORF_li)), rep("overlapping", nrow(overlapping_uORFs_li))))

write.table(all_uORFs_lv, file= paste("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/RI_lv_All_uORFs.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F) 
write.table(all_uORFs_li, file= paste("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/RI_li_All_uORFs.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F) 


all_orfs_bed_lv = read.table("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/RiboTaper/results/pooled_all_lv/translated_ORFs_sorted.bed", sep="\t", header=F, stringsAsFactors=F)
all_orfs_bed_li = read.table("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/RiboTaper/results/pooled_all_liver/translated_ORFs_sorted.bed", sep="\t", header=F, stringsAsFactors=F)
all_orfs_bed_lv$ORF_id_tr = str_split_fixed(all_orfs_bed_lv$V4, ";", 3)[,1]
uorfs_bed_lv = all_orfs_bed_lv[which(as.character(all_orfs_bed_lv$ORF_id_tr) %in% as.character(all_uORFs_lv$ORF_id_tr)),]
all_orfs_bed_li$ORF_id_tr = str_split_fixed(all_orfs_bed_li$V4, ";", 3)[,1]
uorfs_bed_li = all_orfs_bed_li[which(as.character(all_orfs_bed_li$ORF_id_tr) %in% as.character(all_uORFs_li$ORF_id_tr)),]


write.table(uorfs_bed_lv[,c(1:6)], file= "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/RiboTaper/results/pooled_all_lv/uORFs.bed", sep="\t", col.names=F, row.names=F, quote=F) 
write.table(uorfs_bed_li[,c(1:6)], file= "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/RiboTaper/results/pooled_all_liver/uORFs.bed", sep="\t", col.names=F, row.names=F, quote=F) 

save(all_uORFs_lv, all_uORFs_li, file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_tables.RData")

####QUANTIFY

library(rtracklayer)
library(stringr)
cols=c("type","source","gene_id","gene_version","transcript_id","exon_number","gene_name","gene_source","gene_biotype","havana_gene",
			"havana_gene_version","transcript_name","transcript_source","transcript_biotype","exon_id","exon_version")
# load new human GTF 
rn6_anno <- import.gff(file.path("/data/huebner2/DATASETS/GENOMES/Ensembl_Genomes/Ens82/rn6/annotation/Rattus_norvegicus.Rnor_6.0.82.gtf"), colnames=cols)

load("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/uORF_tables.RData")

all_uorfs_bed_lv = data.frame(all_orfs_bed_lv, all_uORFs_lv[match(all_orfs_bed_lv$ORF_id_tr, all_uORFs_lv$ORF_id_tr), c("gene_id", "gene_symbol", "annotation", "uORF_type")])
all_uorfs_bed_li = data.frame(all_orfs_bed_li, all_uORFs_li[match(all_orfs_bed_li$ORF_id_tr, all_uORFs_li$ORF_id_tr), c("gene_id", "gene_symbol", "annotation", "uORF_type")])



uORF_ranges_lv = GRanges(seqnames = all_uorfs_bed_lv$V1, ranges=IRanges(start = all_uorfs_bed_lv$V2, end = all_uorfs_bed_lv$V3), strand = all_uorfs_bed_lv$V6, type = "CDS", source = "Ribotaper", gene_id = paste(all_uorfs_bed_lv$gene_id, all_uorfs_bed_lv$uORF_type, sep="_"), gene_version = 1, transcript_id=all_uorfs_bed_lv$ORF_id_tr, exon_number = NA, gene_name = all_uorfs_bed_lv$gene_symbol, gene_source = NA, gene_biotype=all_uorfs_bed_lv$uORF_type, havanna_gene = NA, havana_gene_version = NA, transcript_name = all_uorfs_bed_lv$ORF_id_tr, transcript_source = NA, transcript_biotype = NA, exon_id = NA, exon_version= NA)

uORF_ranges_li = GRanges(seqnames = all_uorfs_bed_li$V1, ranges=IRanges(start = all_uorfs_bed_li$V2, end = all_uorfs_bed_li$V3), strand = all_uorfs_bed_li$V6, type = "CDS", source = "Ribotaper", gene_id = paste(all_uorfs_bed_li$gene_id, all_uorfs_bed_li$uORF_type, sep="_"), gene_version = 1, transcript_id=all_uorfs_bed_li$ORF_id_tr, exon_number = NA, gene_name = all_uorfs_bed_li$gene_symbol, gene_source = NA, gene_biotype=all_uorfs_bed_li$uORF_type, havanna_gene = NA, havana_gene_version = NA, transcript_name = all_uorfs_bed_li$ORF_id_tr, transcript_source = NA, transcript_biotype = NA, exon_id = NA, exon_version= NA)


###first try quantify uORFs and mORFs seperately using the complete coordinates

#1. all uORFs and mORFs in one gtf file --> non_overlapping part of uORF and mORF, plus all real uORFs
#2. all overlapping uORFs in one gtf --> to get all counts in the overlapping uORFs
#3. CDS counts for all genes - done, take from count matrices

rn6_anno.filt = rn6_anno[which(rn6_anno$type == "CDS"),]

export.gff(rn6_anno.filt, "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/uORFs_quant_complete_lv.gtf")
export.gff(uORF_ranges_lv, "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/uORFs_quant_complete_lv.gtf", append=TRUE)


export.gff(rn6_anno.filt, "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/uORFs_quant_complete_li.gtf")
export.gff(uORF_ranges_li, "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/uORFs_quant_complete_li.gtf", append=TRUE)


export.gff(uORF_ranges_lv[which(uORF_ranges_lv$gene_biotype == "overlapping"),], "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/uORFs_quant_only_overlapping_uORFs_lv.gtf")
export.gff(uORF_ranges_li[which(uORF_ranges_li$gene_biotype == "overlapping"),], "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/uORFs_quant_only_overlapping_uORFs_li.gtf")


####counts

cd /data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/lv_ribo

for i in /data/huebner2/ANALYSES/20160310_svh_SH_riboseq_HxB/results/20160310_T210_b2/*_lv/accepted_hits.bam
  do
    name=$(echo $i | awk -F "/" '{print $8}' )
   qsub -b y -V -cwd -N ${name}_cds -l h_vmem=10G  -l h_rt=05:00:00 "htseq-count -f bam -r pos --stranded=yes -t CDS -i gene_id  \
         $i /data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/uORFs_quant_only_overlapping_uORFs_lv.gtf \
         > /data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/lv_ribo/${name}_CDS_counts.txt"
  done

cd /data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/li_ribo

for i in /data/huebner2/ANALYSES/20160310_svh_SH_riboseq_HxB/results/20160310_T210_b2/*_li/accepted_hits.bam
  do
    name=$(echo $i | awk -F "/" '{print $8}' )
   qsub -b y -V -cwd -N ${name}_cds -l h_vmem=10G  -l h_rt=05:00:00 "htseq-count -f bam -r pos --stranded=yes -t CDS -i gene_id  \
         $i /data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/uORFs_quant_only_overlapping_uORFs_li.gtf \
         > /data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/li_ribo/${name}_CDS_counts.txt"
  done


cd /data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/lv_ribo_indep_uORF

for i in /data/huebner2/ANALYSES/20160310_svh_SH_riboseq_HxB/results/20160310_T210_b2/*_lv/accepted_hits.bam
  do
    name=$(echo $i | awk -F "/" '{print $8}' )
   qsub -b y -V -cwd -N ${name}_cds -l h_vmem=10G  -l h_rt=05:00:00 "htseq-count -f bam -r pos --stranded=yes -t CDS -i gene_id  \
         $i /data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/uORFs_quant_complete_lv.gtf \
         > /data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/lv_ribo_indep_uORF/${name}_CDS_counts.txt"
  done

cd /data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/li_ribo_indep_uORF

for i in /data/huebner2/ANALYSES/20160310_svh_SH_riboseq_HxB/results/20160310_T210_b2/*_li/accepted_hits.bam
  do
    name=$(echo $i | awk -F "/" '{print $8}' )
   qsub -b y -V -cwd -N ${name}_cds -l h_vmem=10G  -l h_rt=05:00:00 "htseq-count -f bam -r pos --stranded=yes -t CDS -i gene_id  \
         $i /data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/uORFs_quant_complete_li.gtf \
         > /data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/Detected_ORFs/quantification_uORFs/li_ribo_indep_uORF/${name}_CDS_counts.txt"
  done



