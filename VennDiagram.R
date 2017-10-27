library(VennDiagram)
library(RColorBrewer)
#display.brewer.all()
fast5<-read.csv(file="Fast5_genelist.txt",header=T)
nudis2<-read.csv(file="Nudis2_genelist.txt",header=T)
teas<-read.csv(file="Teas_genelist.txt",header=T)


geneLS<-c(fast5,nudis2,teas)

#Print last few entries
lapply(geneLS,tail)
#check names
names(geneLS)
#Make VennDiagram with colorbrewer pretty colors
venn.plot<-venn.diagram(geneLS,NULL,fill=brewer.pal(3,"Pastel2"))
grid.draw(venn.plot)
pdf("VennDiagramBaits.pdf")
plot<-grid.draw(venn.plot)
dev.off()





