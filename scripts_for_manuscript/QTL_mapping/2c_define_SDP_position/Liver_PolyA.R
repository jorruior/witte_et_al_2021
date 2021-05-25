#!/usr/bin/Rscript

library(GenomicRanges)
library(eQTLpipeline)
library(data.table)
source("Process_CountData.R")
dir="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Result_tables/"

start =  as.numeric(commandArgs(TRUE)[2])
end =  as.numeric(commandArgs(TRUE)[3])
counter = as.numeric(commandArgs(TRUE)[1])


QTL_cat_test = function(qtl_results, ids_unknown, local_sdps)
{
	zeit2=Sys.time()	
	kind = rep("NA", nrow(qtl_results))
	local_sdps = data.table(local_sdps)
	setkey(local_sdps, chr, start, ID, local_ID)
	local = apply(qtl_results[ids_unknown,], 1, function(x) 
	{
		local_sdps_per_assoc = local_sdps[ID == x[1]]
		local_sdps_per_assoc = data.table(local_sdps_per_assoc)
		setkey(local_sdps_per_assoc, chr)
		new_local = local_sdps_per_assoc[which.min(sapply(local_sdps_per_assoc$local_ID, function(local_id) 
		{		
			ids = c((local_sdps[,.I[local_ID == local_id]]-6):(local_sdps[,.I[local_ID == local_id]]+5))
			ids = ids[ids > 0]
			ids = local_sdps[ids]$ID
			mean(as.numeric(qtl_results[SNP %in% ids & gene == x[2]][[5]])) - as.numeric(x[5])	
		}))]

		return(c(new_local$chr, new_local$local_ID))
		
	})
	zeit4=Sys.time()
	print(zeit4-zeit2)
	
return(local)


}


load(file="/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/QTL_mapping/Result_tables/local_sdps.RData")
load(paste(dir, "20180213_Liver_polyA/image.RData", sep=""))
load(file=paste(dir, "20180213_Liver_polyA/eqtls_sorted.RData", sep=""))
eQTL_results_Liver_polyA  = data.table(eQTL_results_liver_polyA)
setkey(eQTL_results_Liver_polyA ,gene,chrom,snp_pos)



locals_unknowns = QTL_cat_test(eQTL_results_Liver_polyA, ids_unknown[start:end], local_sdps)
save(locals_unknowns, file =  paste(dir, "20180213_Liver_polyA/unknowns_", counter, ".RData", sep=""))

