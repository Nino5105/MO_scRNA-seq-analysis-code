# 0.libaray packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggsci)
library(clusterProfiler)
library(org.Mm.eg.db)
options(future.globals.maxSize = 20000 * 1024^2)

data <- readRDS("Myeloid_anno.rds")
data

metadata <- read.csv("Myeloid_merge_meta.csv",row.names = 1) 
head(data@meta.data)
table(metadata$cell_subtype4)
metadata$cell_subtype5 <- metadata$cell_subtype4

metadata$cell_subtype4 <- factor(metadata$cell_subtype4, levels = rev(c(
  "Mono_C1_Cd14","Mono_C2_Lyz2",
  "Mono_C3_S100a8","Mono_C4_Cd74",
  "Mono_C5_Fcgr4","Megakaryocyte",
  
  "Kupffer_C1_Mrc1","Kupffer_C2_Marco",
  "Kupffer_C3_Itgax","Kupffer_C4_Ccr2",
  "LCM_C1_Marco","LCM_C2_Itgax",
  "LCM_C3_Ccr2","LCM_C4_Mki67",
  
  "SPM_C1_C1qc","SPM_C2_Mrc1","SPM_C3_Marco",
  "SPM_C4_Itgam","SPM_C5_Ccr2","SPM_C6_Ccl5",
  
  "Neutro_Act","Neutro_Home","Platelet","Mast cell","Dendritic cell")))
table(metadata$cell_subtype4 == metadata$cell_subtype5)

data@meta.data <- metadata


DimPlot(data,pt.size = 0.1,group.by = "cell_subtype3")
DimPlot(data,pt.size = 0.1,group.by = "cell_subtype4",label = T)


FeaturePlot(data,features = "Cd68")
FeaturePlot(data,features = "Cd68")


cols <- c("Mono_C1_Cd14" = "#925E9FFF",
          "Mono_C2_Lyz2" = "#0099B4FF",
          "Mono_C3_S100a8" = "#42B540FF",
          "Mono_C4_Cd74" = "#ED0000FF",
          "Mono_C5_Fcgr4" = "#00468BFF",
          "Megakaryocyte" = "#6A6599FF",
          
          "Kupffer_C1_Mrc1" = "#B09C85FF",
          "Kupffer_C2_Marco" = "#7E6148FF",
          "Kupffer_C3_Itgax" = "#DC0000FF",
          "Kupffer_C4_Ccr2" = "#91D1C2FF",
          "LCM_C1_Marco" = "#8491B4FF",
          "LCM_C2_Itgax" = "#F39B7FFF",
          "LCM_C3_Ccr2" = "#3C5488FF",
          "LCM_C4_Mki67" = "#00A087FF",
          
          "SPM_C1_C1qc" = "#4DBBD5FF",
          "SPM_C2_Mrc1" = "#E64B35FF",
          "SPM_C3_Marco" = "#AD002AFF",
          "SPM_C4_Itgam" = "#1B1919FF",
          "SPM_C5_Ccr2" = "#ADB6B6FF",
          "SPM_C6_Ccl5" = "#FDAF91FF",
          
          "Neutro_Act" = "#20854EFF",
          "Neutro_Home" = "#7876B1FF",
          "Platelet" = "#6F99ADFF",
          "Mast cell" = "#0072B5FF",
          "Dendritic cell" = "#BC3C29FF"
          )


DimPlot(data,pt.size = 0.1,group.by = "cell_subtype4",
        label = F,cols = cols,split.by = "orgin")

# subset
data_liver <- subset(data,subset = orgin == "Liver")
data_pbmc <- subset(data,subset = orgin == "PBMC")
data_spleen <- subset(data,subset = orgin == "Spleen")

DimPlot(data_liver,pt.size = 0.1,group.by = "cell_subtype4",
        label = F,cols = cols,split.by = "orgin")
DimPlot(data_pbmc,pt.size = 0.1,group.by = "cell_subtype4",
        label = F,cols = cols,split.by = "orgin")
DimPlot(data_spleen,pt.size = 0.1,group.by = "cell_subtype4",
        label = F,cols = cols,split.by = "orgin")

table(data$orgin,data$cell_subtype4)
table(data$cell_subtype4)

data$orgin <- factor(data$orgin,levels = c("PBMC","Liver","Spleen"))


DoHeatmap(data,group.by = "cell_subtype4",c(
  "Cd68",
  "Cpa3","Cd63",# Mast
  "S100a8","S100a9",# Neutro
  "Ltf","Camp","Ngp",# Neutro-Home
  "Il1b","Csf3r","Ifitm1",# Neutro-act
  "Siglech",# pDC
  "Gp9","Pf4" # Platelet
)) 


gene_pbmc <- c(
  
  "Cd68","Cd14",
  "Lyz2","Cx3cr1",
  "S100a8",
  "Cd74","Fcgr4",
  "Pf4","Ppbp", 
  # 
  # "Il1a","Il1b",
  
  "S100a8","S100a9",# Neutro
  "Csf3r","Ifitm1","Il1b",# Neutro-act
  "Ltf","Camp",#"Ngp",# Neutro-Home
  
  "Pf4","Gp9", # Platelet
  "Cpa3","Cd63","Ccl3",# Mast
  "Siglech" # pDC
)

data_pbmc@active.ident <- data_pbmc$cell_subtype4
p1 <- VlnPlot(data_pbmc, 
              features = gene_pbmc, 
              pt.size=0, fill.by = "ident",
              cols = rev(cols),
              stack = T)  + xlab(NULL) + ylab(NULL) + NoLegend()
p1


# modify_vlnplot <- function(obj,
#                            feature,
#                            pt.size = 0,
#                            plot.margin = unit(c(-0.75, 0, -0.75, 0), "mm"),
#                            ...) {
#   p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
#     xlab("") + ylab(feature) + ggtitle("") +
#     theme(legend.position = "none",
#           axis.text.x = element_blank(),
#           axis.text.y = element_blank(),
#           axis.ticks.x = element_blank(),
#           axis.ticks.y = element_line(),
#           axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
#           plot.margin = plot.margin )
#   return(p)
# }
# 
# StackedVlnPlot<- function(obj, features,
#                           pt.size = 0,
#                           plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
#                           ...) {
#   
#   plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
#   plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
#     theme(axis.text.x=element_text(), axis.ticks.x = element_line())
#   
#   p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
#   return(p)
# }


gene <- c(
  "Cd68","Cd14",
  "Lyz2",#"Cx3cr1",
  "S100a8",#"Ccl3",
  "Cd74","Fcgr4",
  "Pf4",#"Ppbp", 
  
  "C1qa","C1qc",
  "Cd80","Cd86", # M1
  "Mrc1","Cd163", # M2
  "Marco","Siglece",
  "Itgax","Ccr2",
  "Mki67","Cdca3",
  
  "Adgre1","Itgam",
  "Ly6c2", "Ccr2",
  "Siglec1","Ccl5",
  
  "S100a8","S100a9",# Neutro
  "Ltf","Camp","Ngp",# Neutro-Home
  "Il1b","Csf3r","Ifitm1",# Neutro-act
  
  "Gp9","Pf4", # Platelet
  "Cpa3","Cd63",# Mast
  "Siglech" # pDC
)



data@active.ident <- data$cell_subtype4
StackedVlnPlot(data, gene, pt.size=0, cols=cols)

gene_liver <- c(
  "Cd68","Cd14",
  "C1qc",#"C1qa",
  "Cd86",#"Cd80", # M1
  "Mrc1","Cd163", # M2
  "Marco","Siglece",
  "Itgax","Ccr2",
  "Mki67",#"Cdca3",
  
  "S100a8","S100a9",# Neutro
  "Csf3r","Ifitm1",#"Il1b",# Neutro-act
  "Ltf","Camp",#"Ngp",# Neutro-Home
  # "Cd15","Cd10","Cd11b","Cd16",
  # "Cxcr2",#"Cxcr4","Cd62l",
  
  "Pf4",#"Gp9", # Platelet
  "Cpa3",#"Cd63",# Mast
  "Siglech" # pDC
)

data_liver@active.ident <- data_liver$cell_subtype4
StackedVlnPlot(data_liver, gene_liver, pt.size=0, cols=cols)

p2 <- VlnPlot(data_liver, 
              features = gene_liver, 
              pt.size=0, fill.by = "ident",
              cols = rev(cols),
              stack = T)  + xlab(NULL) + ylab(NULL) + NoLegend()
p2





gene_spleen<- c(
  "Cd68","Cd14",
  "Marco","C1qc",
  "Mrc1","Cd163",
  "Itgam","Ly6c2", 
  "Ccr2","Ccl5",
  
  "S100a9",# "S100a8",# Neutro
  "Il1b","Csf3r",#"Ifitm1",# Neutro-act
  "Ltf","Camp",#"Ngp",# Neutro-Home
  
  "Gp9","Pf4", # Platelet
  "Cpa3","Cd63",# Mast
  "Siglech" # pDC
)

data_spleen@active.ident <- data_spleen$cell_subtype4

p3 <- VlnPlot(data_spleen, 
              features =  gene_spleen, 
              pt.size=0, fill.by = "ident",
              cols = rev(cols),
              stack = T) + ylab(NULL) + NoLegend()
p3


p1 / p2 / p3 # 9*13

data@active.ident <- factor(data$cell_subtype4)
markers <- FindAllMarkers(data)
table(markers$cluster)
write.csv(markers,"subtype_markers.csv")

markers2 <- markers

markers2$cluster <- factor(markers2$cluster, levels = c(
  "Mono_C1_Cd14","Mono_C2_Lyz2",
  "Mono_C3_S100a8","Mono_C4_Cd74",
  "Mono_C5_Fcgr4","Megakaryocyte",
  
  "Kupffer_C1_Mrc1","Kupffer_C2_Marco",
  "Kupffer_C3_Itgax","Kupffer_C4_Ccr2",
  "LCM_C1_Marco","LCM_C2_Itgax",
  "LCM_C3_Ccr2","LCM_C4_Mki67",
  
  "SPM_C1_C1qc","SPM_C2_Mrc1","SPM_C3_Marco",
  "SPM_C4_Itgam","SPM_C5_Ccr2","SPM_C6_Ccl5",
  
  "Neutro_Act","Neutro_Home","Platelet","Mast cell","Dendritic cell"))
table(markers2$cluster)

markers_list <- split(markers2$gene,markers2$cluster)

markers_list_GO <- compareCluster(markers_list,
                                  fun="enrichGO",
                                  OrgDb ="org.Mm.eg.db", 
                                  pvalueCutoff=0.01,
                                  keyType ="SYMBOL")


markers_list_GO@compareClusterResult
ggplot(markers_list_GO,aes(x=Cluster,y=Description)) +
  geom_point(aes(size=Count,color=p.adjust)) +
  theme_bw() + scale_color_continuous(low='#E24D36',high='black') + 
  # ggtitle("Up-regulated pathways (Model vs.Control)") +
  scale_size_continuous(name="GeneRatio",range = c(1,10))+
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(colour = "black",size = 12, vjust =1 ),
        axis.text.y = element_text(colour = "black",size = 12, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),color = "black",size = 12),
        axis.title.y = element_text(angle=90)) 

markers_list_GO@compareClusterResult$Cluster <- factor(markers_list_GO@compareClusterResult$Cluster,
                                                       levels = rev(c(
                                                         "Mono_C1_Cd14","Mono_C2_Lyz2",
                                                         "Mono_C3_S100a8","Mono_C4_Cd74",
                                                         "Mono_C5_Fcgr4","Megakaryocyte",
                                                         
                                                         "Kupffer_C1_Mrc1","Kupffer_C2_Marco",
                                                         "Kupffer_C3_Itgax","Kupffer_C4_Ccr2",
                                                         "LCM_C1_Marco","LCM_C2_Itgax",
                                                         "LCM_C3_Ccr2","LCM_C4_Mki67",
                                                         
                                                         "SPM_C1_C1qc","SPM_C2_Mrc1","SPM_C3_Marco",
                                                         "SPM_C4_Itgam","SPM_C5_Ccr2","SPM_C6_Ccl5",
                                                         
                                                         "Neutro_Act","Neutro_Home","Platelet","Mast cell","Dendritic cell")))
ggplot(markers_list_GO,aes(x=Description,y=Cluster)) +
  geom_point(aes(size=Count,color=p.adjust)) +
  theme_bw() + scale_color_continuous(low='#E24D36',high='black') + 
  # ggtitle("Up-regulated pathways (Model vs.Control)") +
  scale_color_gradientn(colours = rev(viridis::viridis(20)), guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"), name = "Average \n expression") +
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(colour = "black",size = 12, vjust =1 ),
        axis.text.y = element_text(colour = "black",size = 12, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),color = "black",size = 12),
        axis.title.y = element_text(angle=90)) 


ggplot(markers_list_GO,aes(x=Cluster,y=Description)) +
  geom_point(aes(size=Count,color=p.adjust)) +
  theme_bw() + scale_color_continuous(low='#E24D36',high='black') + 
  # ggtitle("Up-regulated pathways (Model vs.Control)") +
  scale_color_gradientn(colours = rev(viridis::viridis(20)), guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"), name = "Average \n expression") +
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(colour = "black",size = 12, vjust =1 ),
        axis.text.y = element_text(colour = "black",size = 12, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),color = "black",size = 12),
        axis.title.y = element_text(angle=90))  
#  
# Mono_C1_Cd14 <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# Mono_C2_Lyz2 <- markers2[markers2$cluster == "Mono_C2_Lyz2",]$gene
# Mono_C3_S100a8 <- markers2[markers2$cluster == "Mono_C3_S100a8",]$gene
# Mono_C4_Cd74 <- markers2[markers2$cluster == "Mono_C4_Cd74",]$gene
# Mono_C5_Fcgr4 <- markers2[markers2$cluster == "Mono_C5_Fcgr4",]$gene
# Megakaryocyte <- markers2[markers2$cluster == "Megakaryocyte",]$gene
# 
# Kupffer_C1_Mrc1 <- markers2[markers2$cluster == "Kupffer_C1_Mrc1",]$gene
# Kupffer_C2_Marco <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# Kupffer_C3_Itgax <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# Kupffer_C4_Ccr2 <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# LCM_C1_Marco <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# LCM_C2_Itgax <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# LCM_C3_Ccr2 <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# LCM_C4_Mki67 <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# 
# SPM_C1_C1qc <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# SPM_C2_Mrc1 <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# SPM_C3_Marco <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# SPM_C4_Itgam <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# SPM_C5_Ccr2 <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# SPM_C6_Ccl5 <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# 
# Neutro_Act <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# Neutro_Home <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# Platelet <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# Mast cell <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
# Dendritic cell <- markers2[markers2$cluster == "Mono_C1_Cd14",]$gene
