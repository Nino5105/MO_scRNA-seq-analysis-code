library(tidyverse)
library(patchwork)
library(SCENIC)
library(foreach)
library(pheatmap)
library(limma)
library(Seurat)

setwd("/public/Count/scRNA/multi_organ/Scenic/plot")
rugln<-c("Ets2_extended","Cebpb","Eomes_extended","Irf8_extended","Ezh2","Ezh2_extended","Ybx1_extended","Setbp1_extended","Hmgb2","Hmgb2_extended"
,"Jund_extended","Fos_extended","Fos","Irf5_extended","Stat3_extended","Irf1_extended","Cebpb_extended"
,"Foxp1_extended","Bcl3","Batf_extended","Batf","Irf8","Prdm1","Ets1_extended","Junb","Spi1","Spic_extended","Junb_extended","Klf2_extended"
,"Irf2_extended","Ets1","Irf7","Mef2c","Klf9_extended","Psmd12_extended","Stat1","Il21_extended","Stat3"
,"Spi1_extended","Irf7_extended","Stat2_extended"
,"Stat1_extended","Maf_extended","Hivep2_extended","Hivep2")
rugln_order<-rugln[order(rugln)]

cell<-c("Kuffer_C3_Itgax","Kuffer_C4_Ccr2","LCM_C3_Ccr2","Mono_C5_Fcgr4","Neutro_Act","Platelet","SPM_C3_Itgam","LCM_C2_Itgax")

cell_scenic<-list(NA)

for(i in cell){
  cell_scenic[[i]][["auc"]]<-paste0("/public/Count/scRNA/multi_organ/Scenic/",i,"_exprMat_DEG/SCENIC/int/3.4_regulonAUC.Rds")
  cell_scenic[[i]][["meta"]]<-paste0("/public/Count/scRNA/multi_organ/Scenic/",i,"_metadata.csv")
}

cell_scenic<-cell_scenic[-which(is.na(cell_scenic)==T)]

for (i in cell){
AUCmatrix <- readRDS(cell_scenic[[i]][["auc"]])
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)

RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
AUCmatrix<-t(AUCmatrix)

meta<-read.csv(cell_scenic[[i]][["meta"]])

rownames(meta)<-meta$X
meta<-meta[,-1]

#get averageAUC
cellInfo<-meta
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$type),function(cells) rowMeans(AUCmatrix[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
AUCmatrix_data<-regulonActivity_byCellType_Scaled

#get name
a<-strsplit(rownames(AUCmatrix_data),"_")
#sapply(a,function(x) x[1])
#a_1<-sapply(a,"[",1)
col<-c("c1","c2","c3")
AUC_data <- as.data.frame(matrix(data = NA,nrow=length(rownames(AUCmatrix_data)),ncol=3,dimnames =list(rownames(AUCmatrix_data),col)))

AUC_data$c1<-sapply(strsplit(rownames(AUCmatrix_data),"_"),"[",1)
AUC_data$c2<-sapply(strsplit(rownames(AUCmatrix_data),"_"),"[",2)
AUC_data$c3<-ifelse(AUC_data$c2=="extended", paste0(AUC_data$c1,"_",AUC_data$c2), AUC_data$c1)

#
AUC_data$name<-rownames(AUC_data)
AUC_data<-AUC_data[,c("name","c3")]
#
AUCmatrix_data<-as.data.frame(AUCmatrix_data)
AUCmatrix_data$name<-rownames(AUCmatrix_data)

data<-merge(AUCmatrix_data,AUC_data,by="name",all=T)
#bulid same data
same_name <- intersect(data$c3,rugln)
sub_data <- data[data$c3 %in% same_name,]
rownames(sub_data)<-sub_data$c3
sub_data<-sub_data[,c("Control","Model","ART")]

#build not same
diffname<-setdiff(rugln,same_name)
col<-c("Control","Model","ART")
sub_data_2 <- as.data.frame(matrix(data = NA,nrow=length(diffname),ncol=3,dimnames =list(diffname,col)))

#merge
merge_data<-rbind(sub_data,sub_data_2)

#order
merge_data<-merge_data[rugln_order,]
#color = colorRampPalette(c(“navy”, “white”, “firebrick3”))(102)
p<-pheatmap(merge_data,cluster_cols=F,cluster_row=F,main=paste0(i,"_averageAUC_heatmap"))
ggsave(paste0(i,"_averageAUC_heatmap.pdf"),p,width=5,height=13)
}
print("done!")
