# 1.1 load packages

library(Seurat)
library(GSEABase)
library(Biobase)
library(limma)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(Cairo)
library(plyr)


options(future.globals.maxSize = 60000 * 1024^2)

setwd("/public/Count/scRNA/multi_organ/2.annotation/plot")
scc_integrated<-readRDS("/public/Count/scRNA/multi_organ/2.annotation/multi_organ_filtered.rds")

p<-DimPlot(scc_integrated,label = T,group.by = "cell_type_filter", label.size = 5,raster=FALSE)
p<-p + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15))

pdf("umap_filtered_final.pdf")
print(p)
dev.off()

meta<-scc_integrated@meta.data

summary<-dplyr::select(meta,"type","cell_type_filter")
summary<-as.data.frame(table(summary))
data <- ddply(summary,"type",transform,percent_number=Freq/sum(Freq)*100)

data$cell_type_filter<-factor(data$cell_type_filter,levels=c("RBC","Hepatocyte","Endo","T/NK cell","B cell","Plasma","Platelet","Macro","Neutro","pDC","Mast cell"))

#draw
p<-ggplot(data,aes(type,percent_number,fill=cell_type_filter))+geom_bar(stat="identity",position="stack")
p<-p+theme(axis.title.x = element_text(size=25,family = "myFont",color = "black",face = "bold"))
p<-p+theme_classic()+geom_col()+theme(text=element_text(size=16))+scale_fill_manual(values=c("RBC" = '#58A4C3',"Hepatocyte" = '#E4C755',"Endo" = 'orange',"T/NK cell" = "#00BFC4","B cell" = '#8C549C',
"Plasma" = '#23452F',"Platelet" = '#AB3282',"Macro" = 'red',"Neutro" = '#BD956A',"pDC" = '#E59CC4',"Mast cell" = '#E63863'))
#p<-p+geom_text(size=3.5,aes(label=round(digit,digits = 1)),position = position_stack(0.5))
p<-p+ggtitle("num")+theme(plot.title = element_text(hjust = 0.5))
p
ggsave("type_num.pdf",p,width = 3.5,height = 5)



summary<-dplyr::select(meta,"orgin","cell_type_filter")
summary<-as.data.frame(table(summary))
data <- ddply(summary,"orgin",transform,percent_number=Freq/sum(Freq)*100)

data$cell_type_filter<-factor(data$cell_type_filter,levels=c("RBC","Hepatocyte","Endo","T/NK cell","B cell","Plasma","Platelet","Macro","Neutro","pDC","Mast cell"))

#draw
p<-ggplot(data,aes(orgin,percent_number,fill=cell_type_filter))+geom_bar(stat="identity",position="stack")
p<-p+theme(axis.title.x = element_text(size=25,family = "myFont",color = "black",face = "bold"))
p<-p+theme_classic()+geom_col()+theme(text=element_text(size=16))+scale_fill_manual(values=c("RBC" = '#58A4C3',"Hepatocyte" = '#E4C755',"Endo" = 'orange',"T/NK cell" = "#00BFC4","B cell" = '#8C549C',
"Plasma" = '#23452F',"Platelet" = '#AB3282',"Macro" = 'red',"Neutro" = '#BD956A',"pDC" = '#E59CC4',"Mast cell" = '#E63863'))
#p<-p+geom_text(size=3.5,aes(label=round(digit,digits = 1)),position = position_stack(0.5))
p<-p+ggtitle("num")+theme(plot.title = element_text(hjust = 0.5))
p
ggsave("orgin_num.pdf",p,width = 3.5,height = 5)


print("done!")
