#preseq make meta plots from library complexity

files_lv_polyA = list.files("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/library_complexity/polyA_trimmed/lv", pattern="yield", full.names=T)
files_liver_polyA = list.files("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/library_complexity/polyA_trimmed/liver", pattern="yield", full.names=T)

files_lv_ribo = list.files("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/library_complexity/ribo/lv", pattern="yield", full.names=T)
files_liver_ribo = list.files("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/library_complexity/ribo/liver", pattern="yield", full.names=T)

polyA_lv_libComplexity=NULL
for(file in files_lv_polyA)
{
	
	content = read.table(file, sep="\t", header=T, stringsAsFactors=F) 
curve_file = gsub("yield_estimates", "c_curve", file)	
	content2 = read.table(curve_file, sep="\t", header=T, stringsAsFactors=F) 
	polyA_lv_libComplexity = rbind(polyA_lv_libComplexity, data.frame(basename(file), content[1:1000,], distinct_seq = content2[match(content[1:1000,"TOTAL_READS"], content2$total_reads),]$distinct_reads))
}

ribo_lv_libComplexity=NULL
for(file in files_lv_ribo)
{
	
	content = read.table(file, sep="\t", header=T, stringsAsFactors=F) 
	curve_file = gsub("yield_estimates", "c_curve", file)	
	content2 = read.table(curve_file, sep="\t", header=T, stringsAsFactors=F) 
	ribo_lv_libComplexity = rbind(ribo_lv_libComplexity, data.frame(basename(file), content[1:1000,], distinct_seq = content2[match(content[1:1000,"TOTAL_READS"], content2$total_reads),]$distinct_reads))
}


polyA_liver_libComplexity=NULL
for(file in files_liver_polyA)
{
	
	content = read.table(file, sep="\t", header=T, stringsAsFactors=F) 
	curve_file = gsub("yield_estimates", "c_curve", file)	
	content2 = read.table(curve_file, sep="\t", header=T, stringsAsFactors=F) 
	polyA_liver_libComplexity = rbind(polyA_liver_libComplexity, data.frame(basename(file), content[1:1000,], distinct_seq = content2[match(content[1:1000,"TOTAL_READS"], content2$total_reads),]$distinct_reads))
}

ribo_liver_libComplexity=NULL
for(file in files_liver_ribo)
{
	
	content = read.table(file, sep="\t", header=T, stringsAsFactors=F) 
	curve_file = gsub("yield_estimates", "c_curve", file)	
	content2 = read.table(curve_file, sep="\t", header=T, stringsAsFactors=F) 
	ribo_liver_libComplexity = rbind(ribo_liver_libComplexity, data.frame(basename(file), content[1:1000,], distinct_seq = content2[match(content[1:1000,"TOTAL_READS"], content2$total_reads),]$distinct_reads))
}




library(ggplot2)
pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/library_complexity/curves_estimated_meta.pdf",width=10)
ggplot(data=polyA_lv_libComplexity, aes(x=TOTAL_READS, y = EXPECTED_DISTINCT, col=as.factor(basename.file.)))+geom_line(linetype="dotted")+ylim(c(0, 7e+07))+theme_classic()+ggtitle("Library Complexity PolyA Left-ventricle")+ theme(legend.position="none")

ggplot(data=ribo_lv_libComplexity, aes(x=TOTAL_READS, y = EXPECTED_DISTINCT, col=as.factor(basename.file.)))+geom_line(linetype="dotted")+ylim(c(0, 2e+07))+theme_classic()+ggtitle("Library Complexity Ribo Left-ventricle")+ theme(legend.position="none")

ggplot(data=polyA_liver_libComplexity, aes(x=TOTAL_READS, y = EXPECTED_DISTINCT, col=as.factor(basename.file.)))+geom_line(linetype="dotted")+ylim(c(0, 7e+07))+theme_classic()+ggtitle("Library Complexity PolyA Liver")+ theme(legend.position="none")

ggplot(data=ribo_liver_libComplexity, aes(x=TOTAL_READS, y = EXPECTED_DISTINCT, col=as.factor(basename.file.)))+geom_line(linetype="dotted")+ylim(c(0, 2e+7))+theme_classic()+ggtitle("Library Complexity Ribo Liver")+ theme(legend.position="none")

dev.off()




library(ggplot2)
library(cowplot)
pdf("/data/huebner2/ANALYSES/20171204_RatHeartTranslatome/results/library_complexity/curves_estimated_individual.pdf",width=15, height=25)

p=list()
x=0
for(i in unique(polyA_lv_libComplexity$basename.file.))
{
	x=x+1	
	p[[x]] = ggplot(data=polyA_lv_libComplexity[which(polyA_lv_libComplexity$basename.file. == i),], aes(x=TOTAL_READS, y = EXPECTED_DISTINCT))+geom_ribbon(aes(ymin=LOWER_0.95CI,ymax=UPPER_0.95CI),alpha=0.3)+geom_line(linetype="dotted") + geom_line(linetype="solid", aes(x=TOTAL_READS, y = distinct_seq))+xlim(c(0, 7e+07)) +ylim(c(0, 7e+07))+theme_classic()+ggtitle(paste(gsub("_pol_yield_estimates.txt", "", i), " PolyA lv", sep=""))+ theme(legend.position="none")
	
}

multiplot(plotlist=p, cols = 3)
p=list()
x=0
for(i in unique(ribo_lv_libComplexity$basename.file.))
{
	x=x+1	
	p[[x]] = ggplot(data=ribo_lv_libComplexity[which(ribo_lv_libComplexity$basename.file. == i),], aes(x=TOTAL_READS, y = EXPECTED_DISTINCT))+geom_ribbon(aes(ymin=LOWER_0.95CI,ymax=UPPER_0.95CI),alpha=0.3)+geom_line(linetype="dotted") + geom_line(linetype="solid", aes(x=TOTAL_READS, y = distinct_seq))+xlim(c(0, 2e+07)) +ylim(c(0, 2e+07))+theme_classic()+ggtitle(paste(gsub("_pol_yield_estimates.txt", "", i), " Ribo lv", sep=""))+ theme(legend.position="none")
	
}

multiplot(plotlist=p, cols = 3)


p=list()
x=0
for(i in unique(polyA_liver_libComplexity$basename.file.))
{
	x=x+1	
	p[[x]] = ggplot(data=polyA_liver_libComplexity[which(polyA_liver_libComplexity$basename.file. == i),], aes(x=TOTAL_READS, y = EXPECTED_DISTINCT))+geom_ribbon(aes(ymin=LOWER_0.95CI,ymax=UPPER_0.95CI),alpha=0.3)+geom_line(linetype="dotted") + geom_line(linetype="solid", aes(x=TOTAL_READS, y = distinct_seq))+xlim(c(0, 7e+07)) +ylim(c(0, 7e+07))+theme_classic()+ggtitle(paste(gsub("_pol_yield_estimates.txt", "", i), " PolyA liver", sep=""))+ theme(legend.position="none")
	
}

multiplot(plotlist=p, cols = 3)
p=list()
x=0
for(i in unique(ribo_liver_libComplexity$basename.file.))
{
	x=x+1	
	p[[x]] = ggplot(data=ribo_liver_libComplexity[which(ribo_liver_libComplexity$basename.file. == i),], aes(x=TOTAL_READS, y = EXPECTED_DISTINCT))+geom_ribbon(aes(ymin=LOWER_0.95CI,ymax=UPPER_0.95CI),alpha=0.3)+geom_line(linetype="dotted") + geom_line(linetype="solid", aes(x=TOTAL_READS, y = distinct_seq))+xlim(c(0, 2e+07)) +ylim(c(0, 2e+07))+theme_classic()+ggtitle(paste(gsub("_pol_yield_estimates.txt", "", i), " Ribo liver", sep=""))+ theme(legend.position="none")
	
}

multiplot(plotlist=p, cols = 3)


dev.off()

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}







