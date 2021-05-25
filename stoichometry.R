##Script to plot the stoichometry of specific genes based on the data from congenic rats and the RI panel
library("plyr")
library("ggplot2")

cong = ndds_congenics_unique[,1:10] #File with counts - congenic rats
c1 = c(0,0,0,0,0,2,2,2,2,2) #Genotypes
ri = ndds_ri_unique[,1:29] #File with counts - RI panel
c2 = c(0,2,0,2,2,0,2,0,2,0,0,2,2,2,2,2,0,2,0,2,0,0,0,0,0,0,0,2,2) #Genotypes
cong = cong[rownames(cong) %in% rownames(ri),]
ri = ri[rownames(ri) %in% rownames(cong),]
congri_counts = cbind(cong,ri)
congri_counts = congri_counts[rownames(congri_counts) %in% thin,]
congri = congri_counts/l
congri = congri * tr
c12 = c(c1,c2)
rownames(congri) = thing
congri = t(t(congri) / apply(congri,2,sum))

congri_endo_ribo = melt(congri[,c12 == 0])
congri_endo_ribo["class"] = "Ribo-seq"

cong = ndds_congenics_unique[,11:20]
ri = ndds_ri_unique[,30:58]
congrn_counts = cbind(cong,ri)
congrn_counts = congrn_counts[rownames(congrn_counts) %in% thin,]
congrn = congrn_counts/l
congrn = congrn * tr
colnames(congrn) = colnames(congri)
rownames(congrn) = thing
congrn = t(t(congrn) / apply(congrn,2,sum))

congri_endo_rna = melt(congrn[,c12 == 0])
congri_endo_rna["class"] = "RNA-seq"


aaa = ggplot(rbind(congri_endo_ribo,congri_endo_rna), aes(x = Var2, y = value, fill = grepl("endo", Var2)) ) + geom_bar(stat = "identity") + facet_grid(class ~ Var1) + ylab("Proportion") + xlab("Congenics + RI samples - Endo") + theme_bw() + theme(axis.text.x=element_blank(),legend.position="none") + geom_hline(yintercept = 1/11, linetype="dashed") + geom_hline(yintercept = 7/11, linetype="dashed")


congri_ctrl_ribo = melt(congri[,c12 == 2])
congri_ctrl_ribo["class"] = "Ribo-seq"

congri_ctrl_rna = melt(congrn[,c12 == 2])
congri_ctrl_rna["class"] = "RNA-seq"

congri_merged = rbind(congri_endo_ribo,congri_endo_rna,congri_ctrl_ribo,congri_ctrl_rna)

a = ddply(congri_ctrl_ribo, .(Var1,class), function(x) mean(x$value))
b = ddply(congri_endo_ribo, .(Var1,class), function(x) mean(x$value))
dfribo = data.frame(names = a$Var1, endo = a$V1, ctrl = b$V1)
rownames(dfribo) = dfribo$names
dfribo$names <- NULL

a = ddply(congri_ctrl_rna, .(Var1,class), function(x) mean(x$value))
b = ddply(congri_endo_rna, .(Var1,class), function(x) mean(x$value))
dfrna = data.frame(names = a$Var1, endo = a$V1, ctrl = b$V1)
rownames(dfrna) = dfrna$names
dfrna$names <- NULL

te_counts_endo = congri_counts[,c12 == 0]/congrn_counts[,c12 == 0]
apply(te_counts_endo,1,mean)
te_counts_endom = melt(te_counts_endo)
te_counts_endom["type"] = "endo"

te_counts_ctrl = congri_counts[,c12 == 2]/congrn_counts[,c12 == 2]
apply(te_counts_ctrl,1,mean)
te_counts_ctrlm = melt(te_counts_ctrl)
te_counts_ctrlm["type"] = "ctrl"

te_counts_m = rbind(te_counts_endom,te_counts_ctrlm)

bbb = ggplot(te_counts_m[te_counts_m$Var1 %in% thin,], aes(x = Var1, y = value, fill = type), alpha = 0.5) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = type), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size = 1.2) + theme_bw() + ggtitle("TE thin filament") + scale_x_discrete(labels = thing) + xlab("TE") + ylab("") + geom_hline(yintercept = 1, linetype = "dashed")


congri_counts_m = melt(congri_counts[,c12 == 0])
congrn_counts_m = melt(congrn_counts[,c12 == 0])
congri_counts_m["type"] = "Ribo-seq"
congrn_counts_m["type"] = "RNA-seq"
congr_counts_m = rbind(congri_counts_m, congrn_counts_m)

ccc = ggplot(congr_counts_m[(congr_counts_m$Var1 == "ENSRNOG00000008536")|(congr_counts_m$Var1 == "ENSRNOG00000033734"),], aes(x = Var1, y = value, fill = type)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(colour = type), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size = 1.2) + theme_bw() + ggtitle("Counts subunits") + scale_x_discrete(labels = c("Actc1","Tnnt2")) + xlab("Normalized counts") + ylab("")

d = ggarrange(aaa, ggarrange(bbb, ccc, ncol = 1, nrow = 2))
d