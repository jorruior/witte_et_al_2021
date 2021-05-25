sapply(c("psych","HTSCluster","rtracklayer","dendextend","plyr", "gProfileR","dplyr","gplots","gridExtra",
         "ggplot2","reshape2","ggpubr","ggbeeswarm","devtools","xlsx","stringr","cowplot"),require,character.only = TRUE)

load(file = "/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/DESeq_congenics/FC_plot_riboOnlyGenesNotLocatedChr3.RData")
cong_ribo = read.table("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/deseq_norm/polyA_cong_lv_Ribo_deseqNorm_RNAANDRIBO.txt", sep="\t", header=T, stringsAsFactors=F)
cong_polyA = read.table("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/data/count_matrices/deseq_norm/polyA_cong_lv_RNA_deseqNorm_RNAANDRIBO.txt", sep="\t", header=T, stringsAsFactors=F)

ribo_only_notChr3_go.filt = ribo_only_notChr3_go[which(ribo_only_notChr3_go$domain %in% c("BP", "CC", "MF")),]
#ribo_only_notChr3_go.filt = ribo_only_notChr3_go.filt[-which(ribo_only_notChr3_go.filt$term.id == "KEGG:00000"),]

top20_GO_terms = ribo_only_notChr3_go.filt[order(ribo_only_notChr3_go.filt$p.value), ][1:20,]

cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")
group.colors <- c(BP = cbbPalette[2], CC = cbbPalette[3], MF =cbbPalette[4], keg = cbbPalette[5])
top20_GO_terms$term.name <- factor(top20_GO_terms$term.name, levels = rev(top20_GO_terms$term.name))
pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/DESeq_congenics/GO_top20.pdf")
ggplot(top20_GO_terms, aes(x=term.name, y=-log10(p.value),fill=domain) ) + geom_col() + coord_flip() + xlab("")+scale_fill_manual(values=group.colors)
dev.off()

cong_ribo.means = data.frame(endo = apply(cong_ribo[,1:5], 1, function(x) {mean(x, na.rm=F)}), ctrl = apply(cong_ribo[,6:10], 1, function(x) {mean(x, na.rm=F)}))
cong_polyA.means = data.frame(endo = apply(cong_polyA[,1:5], 1, function(x) {mean(x, na.rm=F)}), ctrl = apply(cong_polyA[,6:10], 1, function(x) {mean(x, na.rm=F)}))


go_term_data_table_FC = data.frame(gene_id = rownames(ribo_only_notChr3), rna_FC = cong_polyA.means[rownames(ribo_only_notChr3),]$ctrl / cong_polyA.means[rownames(ribo_only_notChr3),]$endo, ribo_FC =  cong_ribo.means[rownames(ribo_only_notChr3),]$ctrl / cong_ribo.means[rownames(ribo_only_notChr3),]$endo)


i=4
for(go.term in top20_GO_terms$term.name)
{

	term_genes = unlist(strsplit(top20_GO_terms[which(top20_GO_terms$term.name == go.term), "intersection"], ","))
	go_term_data_table_FC = data.frame(go_term_data_table_FC, term = go_term_data_table_FC$gene_id %in% term_genes)
	colnames(go_term_data_table_FC)[i] = go.term
	i=i+1
}

#go_term_data_table_FC$gene_id = factor(go_term_data_table_FC$gene_id, levels = go_term_data_table_FC[order(log2(go_term_data_table_FC$ribo_FC) - log2(go_term_data_table_FC$rna_FC)),"gene_id"])
go_term_data_table_FC$gene_id = factor(go_term_data_table_FC$gene_id, levels = go_term_data_table_FC[order(log2(go_term_data_table_FC$ribo_FC)),"gene_id"])


pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/DESeq_congenics/FC_GoTerm_plot_commonNorm2.pdf", width=15, height=20)

a = ggplot(go_term_data_table_FC, aes(x=gene_id, y=log2(ribo_FC),fill=log2(ribo_FC))) + 
    geom_bar(stat='identity', width=.7)+ylab("log2(Ribo FC)")+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), legend.key.size = unit(0.3,"cm"),
          legend.title =element_blank(),  legend.text = element_text(size=10),
          plot.margin = unit(c(0,1,0,1), "lines"))+
scale_fill_gradient2(breaks = c(-2.5,0,2.5), low=c("darkblue","blue"), high=c("red","darkred"), 
                         guide="colorbar",limits=c(-2.5,2.5))+ylim(-2.5,2.5)

b = ggplot(go_term_data_table_FC, aes(x=gene_id, y=log2(rna_FC),fill=log2(rna_FC))) + 
    geom_bar(stat='identity', width=.7)+ylab("log2(RNA FC)")+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), legend.key.size = unit(0.3,"cm"),
          legend.title =element_blank(),  legend.text = element_text(size=10),
          plot.margin = unit(c(0,1,0,1), "lines"))+
scale_fill_gradient2(breaks = c(-2.5,0,2.5), low=c("darkblue","blue"), high=c("red","darkred"), 
                         guide="colorbar",limits=c(-2.5,2.5))+ylim(-2.5,2.5)


c = ggplot(go_term_data_table_FC, aes(x=gene_id,y=1, fill=extracellular.matrix))+ geom_bar(stat='identity', width=.7)+ylab("ECM genes")
d = ggplot(go_term_data_table_FC, aes(x=gene_id,y=1, fill=ribosome)) + geom_bar(stat='identity', width=.7)+ylab("Ribosome genes")
e = ggplot(go_term_data_table_FC, aes(x=gene_id,y=1, fill=vesicle)) + geom_bar(stat='identity', width=.7)+ylab("Vesicle genes")


grid.arrange(c,d,e, a, b,  nrow = 5)

dev.off()

save()
