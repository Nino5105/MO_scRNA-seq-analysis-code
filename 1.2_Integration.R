# 1.1 loading packages

library(harmony)
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ggsci)
options(future.globals.maxSize = 60000 * 1024^2)

# 1.2 loading ST-seq datasets

setwd("/public/Count/scRNA/multi_organ/1.cluster/merge_sct")
rds<-readRDS("multi_tissue_scRNA_integrated_Umap.rds")

seurat_obj<-rds
seurat_obj <- seurat_obj %>%
  RunHarmony("sample", assay.use="SCT")

# 1.3 UMAP and clustering with harmonized PCs
seurat_obj <- RunUMAP(seurat_obj, reduction='harmony', dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj, reduction='harmony')
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)

saveRDS(seurat_obj, file="seurat_object_sct_harmony.rds")

p1 <- DimPlot(seurat_obj, group.by='type', reduction='umap') +ggtitle('seurat clusters')
pdf("umap_clusters_sct_harmony_1.pdf", width=7, height=7)
print(p1)
dev.off()

p1 <- DimPlot(seurat_obj, group.by='orgin', reduction='umap') +ggtitle('seurat clusters')
pdf("umap_clusters_sct_harmony_2.pdf", width=7, height=7)
print(p1)
dev.off()

# setwd("/public/Count/scRNA/multi_organ/1.cluster/merge_sct")
# rds<-readRDS("seurat_object_sct_harmony.rds")

# seurat_obj <- RunTSNE(rds, reduction='harmony', dims = 1:30)

# p1 <- DimPlot(seurat_obj, group.by='orgin', reduction='tsne')
# pdf("tsne_clusters_sct_harmony_2.pdf", width=7, height=7)
# print(p1)
# dev.off()

# p1 <- DimPlot(seurat_obj, group.by='type', reduction='tsne')
# pdf("tsne_clusters_sct_harmony_1.pdf", width=7, height=7)
# print(p1)
# dev.off()
