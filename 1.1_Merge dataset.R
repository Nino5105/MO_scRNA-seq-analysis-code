# 1.1 load packages
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ggsci)
options(future.globals.maxSize = 60000 * 1024^2)

# 1.2 load 10X genomics datasets
setwd("/public/Count/scRNA/multi_organ/1.cluster/merge_sct")

samplename = c(
  "/public/QC/PC1/outs/filtered_feature_bc_matrix","/public/QC/PC2/outs/filtered_feature_bc_matrix","/public/QC/PC3/outs/filtered_feature_bc_matrix",
  "/public/QC/PM1/outs/filtered_feature_bc_matrix","/public/QC/PM2/outs/filtered_feature_bc_matrix","/public/QC/PM3/outs/filtered_feature_bc_matrix",
  "/public/QC/PA1/outs/filtered_feature_bc_matrix","/public/QC/PA2/outs/filtered_feature_bc_matrix","/public/QC/PA3/outs/filtered_feature_bc_matrix",

  "/data/rawdata/chenjiayun/02_Chen_AAI_Liver_scRNA/2.Cellranger_results/LC1/outs/filtered_feature_bc_matrix",
  "/data/rawdata/chenjiayun/02_Chen_AAI_Liver_scRNA/2.Cellranger_results/LC2/outs/filtered_feature_bc_matrix",
  "/data/rawdata/chenjiayun/02_Chen_AAI_Liver_scRNA/2.Cellranger_results/LC3/outs/filtered_feature_bc_matrix",
  "/public/QC/LM1/outs/filtered_feature_bc_matrix", "/public/QC/LM2/outs/filtered_feature_bc_matrix","/public/QC/LM3/outs/filtered_feature_bc_matrix",
  "/public/QC/LA1/outs/filtered_feature_bc_matrix", "/public/QC/LA2/outs/filtered_feature_bc_matrix", "/public/QC/LA3/outs/filtered_feature_bc_matrix",

  "/data/rawdata/chenjiayun/10_Chen_MO_scRNA/1.Expression_data/spleen/SC1/outs/filtered_feature_bc_matrix",
  "/data/rawdata/chenjiayun/10_Chen_MO_scRNA/1.Expression_data/spleen/SC2/SC2/outs/filtered_feature_bc_matrix",
  "/data/rawdata/chenjiayun/10_Chen_MO_scRNA/1.Expression_data/spleen/SC3/SC3/outs/filtered_feature_bc_matrix",
  "/public/QC/SM1/outs/filtered_feature_bc_matrix","/public/QC/SM2/outs/filtered_feature_bc_matrix","/public/QC/SM3/outs/filtered_feature_bc_matrix",
  "/public/QC/SA1/outs/filtered_feature_bc_matrix","/public/QC/SA2/outs/filtered_feature_bc_matrix","/public/QC/SA3/outs/filtered_feature_bc_matrix")

names(samplename) = c('PC1','PC2','PC3','PM1','PM2','PM3','PA1','PA2','PA3',
                      'LC1','LC2','LC3','LM1','LM2','LM3','LA1','LA2','LA3',
                      'SC1','SC2','SC3','SM1','SM2','SM3','SA1','SA2','SA3')

print(samplename)

print("1.2 load 10X genomics datasets done!")

# 1.3 read matrix
proj <- list()
for(i in 1:length(samplename)){
  print(names(samplename[i]))
  counts <- Read10X(data.dir = samplename[i])
  newproj <- CreateSeuratObject(counts = counts, min.cells = 3, project = names(samplename[i]))
  #print(names(samplename[i]))
  #print(newproj) # print datasets
  newproj$sample <- names(samplename[i])
  newproj[["percent.mt"]] <- PercentageFeatureSet(object = newproj, pattern = "mt-")
  newproj[["percent.hba"]] <- PercentageFeatureSet(object = newproj, pattern = "Hba-")
  newproj[["percent.hbb"]] <- PercentageFeatureSet(object = newproj, pattern = "Hbb-")
  #saveRDS(newproj,file=paste0(names(samplename[i]),".rds"))
  proj[[names(samplename[i])]] <- newproj
}

#save(proj,file="proj.rdata")
print("1.3 Perform SCTransform on multiple samples done!")

# 1.4 merge
scc_integrated<-merge(proj[[1]],y=c(proj[[2]],proj[[3]],proj[[4]],proj[[5]],proj[[6]],proj[[7]],proj[[8]],proj[[9]],proj[[10]],
                                    proj[[11]],proj[[12]],proj[[13]],proj[[14]],proj[[15]],proj[[16]],proj[[17]],proj[[18]],proj[[19]],
                                    proj[[20]],proj[[21]],proj[[22]],proj[[23]],proj[[24]],proj[[25]],proj[[26]],proj[[27]]))
print(dim(scc_integrated))
print("1.4!")

#1.5 plot data quality before
pdf("multi_tissue_scRNA_before_filter.pdf",height = 5,width = 15)
VlnPlot(scc_integrated, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.hba","percent.hbb"), ncol = 5,pt.size = 0)
dev.off()

#1.6 plot data quality after
scc_integrated <- subset(x = scc_integrated, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 30000 & percent.mt < 15)

pdf("multi_tissue_scRNA_after_filter.pdf",height = 5,width = 15)
VlnPlot(scc_integrated, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.hba","percent.hbb"), ncol = 5,pt.size = 0)
dev.off()

print(dim(scc_integrated))
print("1.6!")

#1.7 sct-pc-umap

scc_integrated<-SCTransform(scc_integrated,return.only.var.genes = FALSE,assay = "RNA",verbose = FALSE)
scc_integrated <- RunPCA(object = scc_integrated,verbose = FALSE)
scc_integrated <- FindNeighbors(scc_integrated,dim=1:20)
scc_integrated <- FindClusters(scc_integrated,resolution = 0.8)
scc_integrated <- RunUMAP(scc_integrated,reduction="pca", dims = 1:20)
print(scc_integrated)

#1.8 add metadata

orgin <- scc_integrated$sample
orgin <- as.character(orgin)
orgin[scc_integrated$sample %in% c('PC1','PC2','PC3','PM1','PM2','PM3','PA1','PA2','PA3')] <- 'PBMC'
orgin[scc_integrated$sample %in% c('LC1','LC2','LC3','LM1','LM2','LM3','LA1','LA2','LA3')] <- 'Liver'
orgin[scc_integrated$sample %in% c('SC1','SC2','SC3','SM1','SM2','SM3','SA1','SA2','SA3')] <- 'Spleen'
table(orgin)

orgin <- factor(orgin,levels = c('PBMC','Liver','Spleen'))
scc_integrated$orgin <- orgin

type <- scc_integrated$sample
type <- as.character(type)
type[scc_integrated$sample %in% c('PC1','PC2','PC3','LC1','LC2','LC3','SC1','SC2','SC3')] <- 'Control'
type[scc_integrated$sample %in% c('PM1','PM2','PM3','LM1','LM2','LM3','SM1','SM2','SM3')] <- 'Model'
type[scc_integrated$sample %in% c('PA1','PA2','PA3','LA1','LA2','LA3','SA1','SA2','SA3')] <- 'ART'
table(type)

type <- factor(type,levels = c("Control","Model","ART"))
scc_integrated$type <- type

#  1.9 saveRDS and plot the Umap
saveRDS(scc_integrated,"multi_tissue_scRNA_integrated_Umap.rds")

multi_tissue_scRNA_integrated_Umap <- DimPlot(scc_integrated,label = TRUE,reduction = "umap")
pdf("multi_tissue_scRNA_integrated_Umap.pdf",height = 10,width = 10)
print(multi_tissue_scRNA_integrated_Umap)
dev.off()

p1 <- DimPlot(scc_integrated,label = TRUE,reduction = "umap",group.by="type")
pdf("multi_tissue_group_Umap.pdf",height = 10,width = 10)
print(p1)
dev.off()

print("done!")

