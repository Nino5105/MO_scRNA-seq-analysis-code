library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(patchwork)
library(cowplot)
library(Cairo)
library(ggplot2)
library(dplyr)

setwd("/public/Count/scRNA/multi_organ/5.interact")

rds<-readRDS("/public/Count/scRNA/multi_organ/2.annotation/multi_organ_filtered.rds")

rds_data<-select(rds@meta.data,cell_type_filter)

rds_data$cell<- unlist(lapply(rds_data$cell_type_filter, function(x) gsub(" ","_",x)))
rds_data$cell<- unlist(lapply(rds_data$cell, function(x) gsub("/","_",x)))

rds_data<-select(rds_data,cell)
rds<-AddMetaData(rds, rds_data, col.name = "cell")

print(table(rds$cell))
print(table(rds$cell_type_filter))

p1 <- DimPlot(rds, reduction = "umap", group.by="cell_type_filter",raster=FALSE)
p2 <- DimPlot(rds, reduction = "umap", group.by="cell", label = TRUE,raster=FALSE,label.size =3.5)
pdf('umap_multi_organ.pdf',width = 10,height = 5)
plot_grid(p1,p2,ncol=2)
dev.off()

Idents(rds)<-rds$cell

meta_data<-rds@meta.data
print(table(meta_data$orgin,meta_data$type))
meta<-list(NA)
data<-list(NA)
for (j in unique(meta_data$orgin)){
        for (i in unique(meta_data$type)){
                meta_sub<-meta_data[which(meta_data$type==i & meta_data$orgin==j),]
                cluster<-subset(rds,cells=rownames(meta_sub))
                p<-DimPlot(cluster,label = FALSE,label.size = 1,cols = c("RBC" = '#58A4C3',"Hepatocyte" = '#E4C755',"Endo" = 'orange',"T/NK cell" = "#00BFC4","B cell" = '#8C549C',
"Plasma" = '#23452F',"Platelet" = '#AB3282',"Macro" = 'red',"Neutro" = '#BD956A',"pDC" = '#E59CC4',"Mast cell" = '#E63863'),group.by = "cell_type_filter",raster=FALSE) + ggtitle(paste0(i,"_",j,"_cluster_umap"))
                p<-CairoTIFF(paste0(i,j,"_cluster_umap.tiff"),units="in",res=200,height=4,width=6)
                print(p)
                dev.off()
                title<-paste0(i,"_",j)
                print(title)
                print(dim(cluster))
                meta[[title]]<-meta_sub
                #data[[title]]<-as.matrix(rds@assays$SCT@data[,as.character(row.names(meta_sub))])
                data[[title]]<-as.matrix(cluster@assays$SCT@data)

}
}
meta<-meta[-which(is.na(meta)==T)]
data<-data[-which(is.na(data)==T)]

#a=c（）；例如变量是i，a=c（a，i）

print(names(meta))
print(names(data))

for (i in names(data)){
        cellchat <- createCellChat(object = data[[i]], meta = meta[[i]], group.by = "cell")
        CellChatDB <- CellChatDB.mouse
        CellChatDB.use <- CellChatDB # simply use the default CellChatDB
        cellchat@DB <- CellChatDB.use
        cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
        future::plan("multiprocess", workers = 4) # do parallel
        cellchat <- identifyOverExpressedGenes(cellchat)
        cellchat <- identifyOverExpressedInteractions(cellchat)
        # project gene expression data onto PPI network (optional)
        cellchat <- projectData(cellchat, PPI.mouse)
        cellchat <- computeCommunProb(cellchat)
        # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
        cellchat <- filterCommunication(cellchat, min.cells = 10)
        #Infer the cell-cell communication at a signaling pathway level
        cellchat <- computeCommunProbPathway(cellchat)
        #Calculate the aggregated cell-cell communication network
        cellchat <- aggregateNet(cellchat)
        groupSize <- as.numeric(table(cellchat@idents))
        pdf(paste0(names(data[i]),"_cellchat_",names(meta[i]),"_interactions.pdf"))
        netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
        dev.off()
        pdf(paste0(names(data[i]),"_cellchat_",names(meta[i]),"_weights.pdf"))
        netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
        dev.off()
        cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
        saveRDS(cellchat, file=paste0(names(data[i]),"_cellchat_",names(meta[i]),".rds"))
}
print("done!")
