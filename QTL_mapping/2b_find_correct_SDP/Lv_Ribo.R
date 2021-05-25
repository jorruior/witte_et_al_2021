library(GenomicRanges)
library(eQTLpipeline)
library(data.table)
source("~/work/RI_translational_efficiency/scripts/R/eQTLpipeline_functions/SeqQTL_oldCluster.R")
dir="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Result_tables/"

#genotype_info = read.table("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Input/snp_pos.txt", sep="\t", header=T, stringsAsFactors=F) 
#genotypes = read.table("~/work/RI_translational_efficiency/data/20160318/RN6_SDPs/HXBBXH-SDP_RN60_LDB-genome.txt",  sep="\t", comment.char = "#", header=F, stringsAsFactors=F, colClasses="character", check.names=F)
#colnames(genotypes) = c("ID", "RI_geno", "SDP_local")
#genotypes = genotypes[match(genotype_info$snp_id, genotypes$ID),]
#local_sdps = read.table("/home/fkreuch/work/RI_translational_efficiency/data/20160318/RN6_SDPs/HXBBXH-SDP_RN60_LDB-local.txt", sep="\t", comment.char = "#", header=F, stringsAsFactors=F, colClasses="character", check.names=F)
#colnames(local_sdps) = c("local_ID", "chr", "start", "end", "Marker_RN60_CSV", "Marker_count")
#rownames(local_sdps) = local_sdps$local_ID
#local_sdps = local_sdps[,c(1:4)]
#SDP_index = apply(local_sdps, 1, function(x) grep(x[1], genotypes$SDP_local))
#SDP_index = SDP_index[sapply(SDP_index, function(x) length(x)!=0)]
#local_sdps = data.frame(local_sdps[match(names(SDP_index), local_sdps$local_ID),], genotypes[as.numeric(SDP_index), ])


#genes_names = gtf2GR("/data/huebner2/DATASETS/GENOMES/Ensembl_Genomes/Ens82/rn6/annotation/#Rattus_norvegicus.Rnor_6.0.82.gtf", c("gene_id", "gene_name", "gene_biotype"))
#genes_names = genes_names[values(genes_names)[,"type"] == "gene"]
#genes_names = genes_names[which(seqnames(genes_names) %in% c(1:20, "X")),]

#local_SDPs_ranges = GRanges(seqnames = local_sdps$chr, ranges=IRanges(start=as.numeric(local_sdps$start), end = as.numeric(local_sdps$end)), global_ID = local_sdps$ID, local_ID = local_sdps$local_ID, add_local_IDs = local_sdps$SDP_local)
#overlap_sdp_gene = findOverlaps(genes_names, local_SDPs_ranges, type="within")
#genes_names$SDP_local = "noSDP"
#genes_names[queryHits(overlap_sdp_gene)]$SDP_local = local_SDPs_ranges[subjectHits(overlap_sdp_gene)]$local_ID
#genes_names$addSDP_local = "noSDP"
#genes_names[queryHits(overlap_sdp_gene)]$addSDP_local = local_SDPs_ranges[subjectHits(overlap_sdp_gene)]$add_local_IDs
#genes_names$SDP_global = "noSDP"
#genes_names[queryHits(overlap_sdp_gene)]$SDP_global = local_SDPs_ranges[subjectHits(overlap_sdp_gene)]$global_ID

#save(local_sdps, file=paste(dir, "local_sdps.RData", sep=""))

QTL_cat_easy = function(qtl_results, local_sdps)
{
	zeit2=Sys.time()	
	kind = rep("NA", nrow(qtl_results))
	library(parallel)
	library(snow)
	#cl= makeCluster(8)
	local = apply(qtl_results, 1, function(x) 
	{
		library(data.table)		
		load(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Result_tables/local_sdps.RData")
		local_sdps = data.table(local_sdps)		
		setkey(local_sdps, ID)		
		local_sdps_per_assoc = local_sdps[ID == x[1]]
		local_sdps_per_assoc = data.table(local_sdps_per_assoc)
		setkey(local_sdps_per_assoc, chr)
		#get information for surrounding SDPs - assume surrounding SNPs of true association is more similar
		if((length(unique(local_sdps_per_assoc$chr)) == 1) | (x[1] == x[14]))
		{
			new_local = local_sdps_per_assoc[1,]
		}
		else
		{
			new_local = data.frame(chr=-1, local_id="xx")
		}
		write(paste(x[2], new_local$chr, new_local$local_ID, sep="\t"), file = paste(dir, "20180213_Lv_Ribo/intermediate_sub_preprocess.txt", sep=""), append=TRUE)
		return(c(new_local$chr, new_local$local_ID))
		
	})
	zeit4=Sys.time()
	print(zeit4-zeit2)
return(local)


}


#load(file=paste(dir, "20180213_Lv_Ribo/eqtls_sorted.RData", sep=""))
#eQTL_results_Lv_Ribo  = data.table(eQTL_results_lv_ribo)
#setkey(eQTL_results_Lv_Ribo ,gene,chrom,snp_pos)


#locals = QTL_cat_easy(eQTL_results_Lv_Ribo , local_sdps)

#save(locals, file =  paste(dir, "20180213_Lv_Ribo/intermediate_prepro.RData", sep=""))
load(file =  paste(dir, "20180213_Lv_Ribo/intermediate_prepro.RData", sep=""))
test = t(as.data.frame(locals))
ids_unknown = which(test[,1] == "-1")


save.image(paste(dir, "20180213_Lv_Ribo/image.RData", sep=""))

