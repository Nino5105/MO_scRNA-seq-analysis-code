# 1.1 load packages
library(Seurat)
library(GSEABase)
library(Biobase)
library(limma)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
options(future.globals.maxSize = 60000 * 1024^2)

# 1.2 load 10X genomics datasets
scc_integrated <- RunPCA(object = scc_integrated,verbose = FALSE)
scc_integrated <- RunUMAP(scc_integrated,reduction="pca", dims = 1:20)
scc_integrated <- FindNeighbors(scc_integrated,dim=1:20)
scc_integrated <- FindClusters(scc_integrated,resolution = 0.8)
scc_integrated <- FindClusters(scc_integrated,resolution = 1.0)

scc_integrated$cell_type <- factor(scc_integrated$cell_type,levels=c("RBC","Hepatocyte","Endo","T/NK cell","B cell","Plasma","Platelet","Macro","Neutro","pDC","Mast cell"))

# 1.3 RenameIdents
scc_integrated <- RenameIdents(scc_integrated,"0"="B cell", "1" = "Macro", "2" = "B cell", "3" = "Endo", "4" = "B cell", "5" = "RBC",
"6" = "T/NK cell", "7" = "RBC", "8" = "Macro","9" = "RBC","10" = "T/NK cell",
"11" = "Macro", "12" = "T/NK cell", "13" = "B cell","14" = "Neutro",
"15" = "RBC", "16" = "RBC", "17" = "B cell", "18" = "T/NK cell", "19" = "T/NK cell",
"20" = "T/NK cell", "21" = "Platelet", "22" = "Endo", "23" = "T/NK cell", "24" = "Plasma",
"25" = "Endo", "26" = "RBC", "27" = "Macro", "28" = "Macro",
"29" = "T/NK cell","30" = "Endo","31" = "Plasma","32" = "Platelet","33" = "B cell",
"34" = "RBC","35" = "Macro","36" = "Neutro","37" = "T/NK cell","38" = "RBC",
"39" = "pDC","40" = "RBC","41" = "T/NK cell","42" = "RBC",
"43" = "Hepatocyte","44" = "Macro","45" = "B cell","46" = "Endo","47" = "Macro",
"48" = "Mast cell","49" = "RBC","50" = "Macro","51" = "RBC","52" = "Endo","53" = "Hepatocyte")

scc_integrated$cell_type_filter<-Idents(scc_integrated)
scc_integrated$cell_type_filter<-factor(scc_integrated$cell_type_filter,levels=c("RBC","Hepatocyte","Endo","T/NK cell","B cell","Plasma","Platelet","Macro","Neutro","pDC","Mast cell"))

# 1.3 visualization

pdf("umap_filtered_final.pdf")
print(p)
dev.off()
p<-DimPlot(scc_integrated,label = T,label.size = 5,cols = c("RBC" = '#58A4C3',"Hepatocyte" = '#E4C755',"Endo" = 'orange',"T/NK cell" = "#00BFC4","B cell" = '#8C549C',
"Plasma" = '#23452F',"Platelet" = '#AB3282',"Macro" = 'red',
"Neutro" = '#BD956A',"pDC" = '#E59CC4',"Mast cell" = '#E63863'),group.by = "cell_type_filter",raster=FALSE)
p<-p + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15))

pdf("umap_filtered_final_2.pdf")
print(p)
dev.off()

Idents(scc_integrated)<-scc_integrated$SCT_snn_res.1
markers<-FindMarkers(scc_integrated,ident.1 = 37, logfc.threshold = 0.25)

p<-DimPlot(scc_integrated,label.size = 5,cols = c('#53A85F','#F1BB72','#E95C59'),group.by = "orgin",raster=FALSE)
p<-p + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15))
pdf("umap_filtered_orgin_final.pdf")
print(p)
dev.off()

p<-DimPlot(scc_integrated,label.size = 5,cols = c("#3B4992","#EE0000","#008B45"),group.by = "type",raster=FALSE)
p<-p + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15))
pdf("umap_filtered_type_final.pdf")
print(p)
dev.off()

saveRDS(scc_integrated,"multi_organ_filtered.rds")



