library(eQTLpipeline)
library(GenomicRanges)

source("Process_CountData.R") ##gtf2GR(), get.htseq.count.matrix()

#load Ensembl GTF to GenomicRanges object
genesGTF = gtf2GR("Rattus_norvegicus.Rnor_6.0.82.gtf", c("gene_id", "gene_name", "gene_biotype"))

#GenomicRanges with gene-centric to derive gene symbols
genes_names = genesGTF[values(genesGTF)[,"type"] == "gene"]
genes_names = genes_names[which(seqnames(genes_names) %in% c(1:20, "X")),]

#GenomicRanges with exon information for non-coding genes (length normalization)
genes = genesGTF[values(genesGTF)[,"type"] == "exon"]
genes = genes[which(seqnames(genes) %in% c(1:20, "X")),]

#GenomicRanges with CDS information for coding genes (length normalization)
genes_CDS = genesGTF[values(genesGTF)[,"type"] == "CDS"]
genes_CDS = genes_CDS[which(seqnames(genes_CDS) %in% c(1:20, "X")),]

#determine all genes that are assigned as protein coding
list_pc = unique(genes[which(genes$gene_biotype == "protein_coding"), ]$gene_id)

#function to make (length normalized) count-matrices 
count_matrices = function(path_full, pattern_gene, pattern_CDS, prefix, list_pc, genes, genes_CDS)
{
	
	
	gene_counts = get.htseq.count.matrix(dir = path_full, pattern=pattern_gene)
	cds_counts  = get.htseq.count.matrix(dir = path_full, pattern= pattern_CDS)
	
	proteincoding_lincRNA = rbind(cds_counts[list_pc, ], gene_counts)
	proteincoding_lincRNA = proteincoding_lincRNA[match(genes_names$gene_id, rownames(proteincoding_lincRNA)),]

	genes_merged = c(genes[which(!(genes$gene_id %in% list_pc)),], genes_CDS[which(genes_CDS$gene_id %in% list_pc),])
	fpkm = get.rpkm(genes_merged, proteincoding_lincRNA)

	return(fpkm)

}

RI_lv_polyA_fpkm = count_matrices("RatHeartTranslatome/data/RI_panel/counts/polyA_trimmed","exon","CDS", "polyA_RI_lv", list_pc, genes, genes_CDS)
RI_lv_ribo_fpkm = count_matrices("RatHeartTranslatome/data/RI_panel/counts/ribo","exon","CDS", "ribo_RI_lv", list_pc, genes, genes_CDS)

RI_li_polyA_fpkm = count_matrices("RatHeartTranslatome/data/RI_panel_liver/polyA_trimmed","exon","CDS", "polyA_RI_li", list_pc, genes, genes_CDS)
RI_li_ribo_fpkm = count_matrices("RatHeartTranslatome/data/RI_panel_liver/ribo","exon","CDS", "ribo_RI_li", list_pc, genes, genes_CDS)

congenics_lv_polyA_fpkm = count_matrices("RatHeartTranslatome/data/congenics/counts/polyA_trimmed","exon","CDS", "polyA_cong_lv", list_pc, genes, genes_CDS)
congenics_lv_ribo_fpkm = count_matrices("RatHeartTranslatome/data/congenics/counts/ribo","exon","CDS", "ribo_cong_lv", list_pc, genes, genes_CDS)

