options(stringsAsFactors = F)

library(seurat)
library(ggradar) # devtools::install_github("ricardo-bion/ggradar",force = TRUE)
library(patchwork)
library(ggplot2)

T_NK <- readRDS("/public/Count/scRNA/multi_organ/2.annotation/subset_dir/T_NK_subset.rds")
T_NK # 42639 features across 48666 samples within 2 assays

Idents(T_NK) <- T_NK$type
DimPlot(T_NK)

T_NK_PBMC <- subset(T_NK,subset = orgin == "PBMC")
T_NK_PBMC@active.ident <- T_NK_PBMC$type
DimPlot(T_NK_PBMC)
T_NK_PBMC_MVC <- T_NK_PBMC_MVC <- FindMarkers(T_NK_PBMC,ident.1 = "Model",ident.2 = "Control",only.pos = F,logfc.threshold = 0)
T_NK_PBMC_MVC$type = ifelse(T_NK_PBMC_MVC$p_val_adj < 0.05 & abs(T_NK_PBMC_MVC$avg_log2FC) >= 0.25, 
                                  ifelse(T_NK_PBMC_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(T_NK_PBMC_MVC$type)
 Down Stable     Up
   400   4629    164

T_NK_PBMC_MVC$orgin <- "PBMC"
write.csv(T_NK_PBMC_MVC,"T_NK_PBMC_MVC.csv")

T_NK_PBMC_MVA <- FindMarkers(T_NK_PBMC,ident.1 = "Model",ident.2 = "ART",only.pos = F,logfc.threshold = 0)
T_NK_PBMC_MVA$type = ifelse(T_NK_PBMC_MVA$p_val_adj < 0.05 & abs(T_NK_PBMC_MVA$avg_log2FC) >= 0.25, 
                                  ifelse(T_NK_PBMC_MVA$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(T_NK_PBMC_MVA$type)
  Down Stable     Up
    96   3859     72

T_NK_PBMC_MVA$orgin <- "PBMC"
write.csv(T_NK_PBMC_MVA,"T_NK_PBMC_MVA.csv")

T_NK_Liver <- subset(T_NK,subset = orgin == "Liver")
T_NK_Liver@active.ident <- T_NK_Liver$type
DimPlot(T_NK_Liver)
T_NK_Liver_MVC <- T_NK_Liver_MVC <- FindMarkers(T_NK_Liver,ident.1 = "Model",ident.2 = "Control",only.pos = F,logfc.threshold = 0)
T_NK_Liver_MVC$type = ifelse(T_NK_Liver_MVC$p_val_adj < 0.05 & abs(T_NK_Liver_MVC$avg_log2FC) >= 0.25, 
                                  ifelse(T_NK_Liver_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(T_NK_Liver_MVC$type)
  Down Stable     Up
   296   5127    833

T_NK_Liver_MVC$orgin <- "Liver"
write.csv(T_NK_Liver_MVC,"T_NK_Liver_MVC.csv")

T_NK_Liver_MVA <- FindMarkers(T_NK_Liver,ident.1 = "Model",ident.2 = "ART",only.pos = F,logfc.threshold = 0)
T_NK_Liver_MVA$type = ifelse(T_NK_Liver_MVA$p_val_adj < 0.05 & abs(T_NK_Liver_MVA$avg_log2FC) >= 0.25, 
                                  ifelse(T_NK_Liver_MVA$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(T_NK_Liver_MVA$type)
Down Stable     Up
    56   6045    129

T_NK_Liver_MVA$orgin <- "Liver"
write.csv(T_NK_Liver_MVA,"T_NK_Liver_MVA.csv")

T_NK_Spleen <- subset(T_NK,subset = orgin == "Spleen")
T_NK_Spleen@active.ident <- T_NK_Spleen$type
DimPlot(T_NK_Spleen)
T_NK_Spleen_MVC <- T_NK_Spleen_MVC <- FindMarkers(T_NK_Spleen,ident.1 = "Model",ident.2 = "Control",only.pos = F,logfc.threshold = 0)
T_NK_Spleen_MVC$type = ifelse(T_NK_Spleen_MVC$p_val_adj < 0.05 & abs(T_NK_Spleen_MVC$avg_log2FC) >= 0.25, 
                                  ifelse(T_NK_Spleen_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(T_NK_Spleen_MVC$type)
  Down Stable     Up
   160   5869    327

T_NK_Spleen_MVC$orgin <- "Spleen"
write.csv(T_NK_Spleen_MVC,"T_NK_Spleen_MVC.csv")

T_NK_Spleen_MVA <- FindMarkers(T_NK_Spleen,ident.1 = "Model",ident.2 = "ART",only.pos = F,logfc.threshold = 0)
T_NK_Spleen_MVA$type = ifelse(T_NK_Spleen_MVA$p_val_adj < 0.05 & abs(T_NK_Spleen_MVA$avg_log2FC) >= 0.25, 
                                  ifelse(T_NK_Spleen_MVA$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(T_NK_Spleen_MVA$type)
  Down Stable     Up
    55   5894     59

T_NK_Spleen_MVA$orgin <- "Spleen"
write.csv(T_NK_Spleen_MVA,"T_NK_Spleen_MVA.csv")

rm(list = ls())


Macro <- readRDS("/public/Count/scRNA/multi_organ/2.annotation/subset_dir/Macro_subset.rds")
Macro # 42639 features across 40902 samples within 2 assays

Idents(Macro) <- Macro$type
DimPlot(Macro)

Macro_PBMC <- subset(Macro,subset = orgin == "PBMC")
Macro_PBMC@active.ident <- Macro_PBMC$type
DimPlot(Macro_PBMC)
Macro_PBMC_MVC <- Macro_PBMC_MVC <- FindMarkers(Macro_PBMC,ident.1 = "Model",ident.2 = "Control",only.pos = F,logfc.threshold = 0)
Macro_PBMC_MVC$type = ifelse(Macro_PBMC_MVC$p_val_adj < 0.05 & abs(Macro_PBMC_MVC$avg_log2FC) >= 0.25, 
                                  ifelse(Macro_PBMC_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Macro_PBMC_MVC$type)
  Down Stable     Up
   423   5618    363

Macro_PBMC_MVC$orgin <- "PBMC"
write.csv(Macro_PBMC_MVC,"Macro_PBMC_MVC.csv")

Macro_PBMC_MVA <- FindMarkers(Macro_PBMC,ident.1 = "Model",ident.2 = "ART",only.pos = F,logfc.threshold = 0)
Macro_PBMC_MVA$type = ifelse(Macro_PBMC_MVA$p_val_adj < 0.05 & abs(Macro_PBMC_MVA$avg_log2FC) >= 0.25, 
                                  ifelse(Macro_PBMC_MVA$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Macro_PBMC_MVA$type)
  Down Stable     Up
   201   4965    212

Macro_PBMC_MVA$orgin <- "PBMC"
write.csv(Macro_PBMC_MVA,"Macro_PBMC_MVA.csv")

Macro_Liver <- subset(Macro,subset = orgin == "Liver")
Macro_Liver@active.ident <- Macro_Liver$type
DimPlot(Macro_Liver)
Macro_Liver_MVC <- Macro_Liver_MVC <- FindMarkers(Macro_Liver,ident.1 = "Model",ident.2 = "Control",only.pos = F,logfc.threshold = 0)
Macro_Liver_MVC$type = ifelse(Macro_Liver_MVC$p_val_adj < 0.05 & abs(Macro_Liver_MVC$avg_log2FC) >= 0.25, 
                                  ifelse(Macro_Liver_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Macro_Liver_MVC$type)
  Down Stable     Up
   348   5042    809

Macro_Liver_MVC$orgin <- "Liver"
write.csv(Macro_Liver_MVC,"Macro_Liver_MVC.csv")

Macro_Liver_MVA <- FindMarkers(Macro_Liver,ident.1 = "Model",ident.2 = "ART",only.pos = F,logfc.threshold = 0)
Macro_Liver_MVA$type = ifelse(Macro_Liver_MVA$p_val_adj < 0.05 & abs(Macro_Liver_MVA$avg_log2FC) >= 0.25, 
                                  ifelse(Macro_Liver_MVA$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Macro_Liver_MVA$type)
  Down Stable     Up
   133   5850    172

Macro_Liver_MVA$orgin <- "Liver"
write.csv(Macro_Liver_MVA,"Macro_Liver_MVA.csv")

Macro_Spleen <- subset(Macro,subset = orgin == "Spleen")
Macro_Spleen@active.ident <- Macro_Spleen$type
DimPlot(Macro_Spleen)
Macro_Spleen_MVC <- Macro_Spleen_MVC <- FindMarkers(Macro_Spleen,ident.1 = "Model",ident.2 = "Control",only.pos = F,logfc.threshold = 0)
Macro_Spleen_MVC$type = ifelse(Macro_Spleen_MVC$p_val_adj < 0.05 & abs(Macro_Spleen_MVC$avg_log2FC) >= 0.25, 
                                  ifelse(Macro_Spleen_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Macro_Spleen_MVC$type)
  Down Stable     Up
   166   4534    603

Macro_Spleen_MVC$orgin <- "Spleen"
write.csv(Macro_Spleen_MVC,"Macro_Spleen_MVC.csv")

Macro_Spleen_MVA <- FindMarkers(Macro_Spleen,ident.1 = "Model",ident.2 = "ART",only.pos = F,logfc.threshold = 0)
Macro_Spleen_MVA$type = ifelse(Macro_Spleen_MVA$p_val_adj < 0.05 & abs(Macro_Spleen_MVA$avg_log2FC) >= 0.25, 
                                  ifelse(Macro_Spleen_MVA$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Macro_Spleen_MVA$type)
Down Stable     Up
   225   5555    103

Macro_Spleen_MVA$orgin <- "Spleen"
write.csv(Macro_Spleen_MVA,"Macro_Spleen_MVA.csv")

rm(list = ls())

#  DEGs analysis

data <- read.csv("DEG_num.csv", check.names=FALSE)
data$Orgin_type <- factor(data$Orgin_type,levels = c("PBMC","Liver","Spleen"))

data_MVC_up <- data[data$DEG_type == "MVC_Up",]

data_MVC_up_P <- ggradar(
  data_MVC_up[,c(2:10)], 
  values.radar = c("0", "500", "1000"),
  grid.min = 0, 
  grid.mid = 500, 
  grid.max = 1000, 
  
  group.line.width = 1.5, 
  group.point.size = 3,
  group.colours = c("#53A85F", "#F1BB72", "#E95C59"),
  
  # background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
) + ggtitle("Up-regulated DEGs numbers \n (Model vs. Control)")
data_MVC_up_P

data_MVC_Down <- data[data$DEG_type == "MVC_Down",]

data_MVC_Down_P <- ggradar(
  data_MVC_Down[,c(2:10)], 
  values.radar = c("0", "500", "1000"),
  grid.min = 0, 
  grid.mid = 500, 
  grid.max = 1000, 
  
  group.line.width = 1.5, 
  group.point.size = 3,
  group.colours = c("#53A85F", "#F1BB72", "#E95C59"),
  
  # background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
) + ggtitle("Down-regulated DEGs numbers \n (Model vs. Control)")

data_MVA_up <- data[data$DEG_type == "MVA_Up",]

data_MVA_up_P <- ggradar(
  data_MVA_up[,c(2:10)], 
  values.radar = c("0", "250", "500"),
  grid.min = 0, 
  grid.mid = 250, 
  grid.max = 500, 
  
  group.line.width = 1.5, 
  group.point.size = 3,
  group.colours = c("#53A85F", "#F1BB72", "#E95C59"),
  
  # background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
) + ggtitle("Up-regulated DEGs numbers \n (Model vs. ART)")
data_MVA_up_P

data_MVA_Down <- data[data$DEG_type == "MVA_Down",]

data_MVA_Down_P <- ggradar(
  data_MVA_Down[,c(2:10)], 
  values.radar = c("0", "250", "500"),
  grid.min = 0, 
  grid.mid = 250, 
  grid.max = 500, 
  
  group.line.width = 1.5, 
  group.point.size = 3,
  group.colours = c("#53A85F", "#F1BB72", "#E95C59"),
  
  # background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
) + ggtitle("Down-regulated DEGs numbers \n (Model vs. ART)")

(data_MVC_up_P + data_MVC_Down_P)/(data_MVA_up_P + data_MVA_Down_P)


# Go enrichment

# T NK cell
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(ggsci)
library(stringr)

T_NK_PBMC_MVC <- read.csv("T_NK/T_NK_PBMC_MVC.csv",row.names = 1,check.names = F)
T_NK_PBMC_MVC_up <- row.names(T_NK_PBMC_MVC[T_NK_PBMC_MVC$type == "Up",])
# T_NK_PBMC_MVC_down <- row.names(T_NK_PBMC_MVC[T_NK_PBMC_MVC$type == "Down",])

T_NK_Liver_MVC <- read.csv("T_NK/T_NK_Liver_MVC.csv",row.names = 1,check.names = F)
T_NK_Liver_MVC_up <- row.names(T_NK_Liver_MVC[T_NK_Liver_MVC$type == "Up",])
# T_NK_Liver_MVC_down <- row.names(T_NK_Liver_MVC[T_NK_Liver_MVC$type == "Down",])

T_NK_Spleen_MVC <- read.csv("T_NK/T_NK_Spleen_MVC.csv",row.names = 1,check.names = F)
T_NK_Spleen_MVC_up <- row.names(T_NK_Spleen_MVC[T_NK_Spleen_MVC$type == "Up",])
# T_NK_Spleen_MVC_down <- row.names(T_NK_Spleen_MVC[T_NK_Spleen_MVC$type == "Down",])

ALL_DEG_up <- list(
  T_NK_PBMC = T_NK_PBMC_MVC_up,
  T_NK_Liver = T_NK_Liver_MVC_up,
  T_NK_Spleen = T_NK_Spleen_MVC_up)

T_NK_MVC_up_GO <- compareCluster(ALL_DEG_up,
                                fun="enrichGO",
                                OrgDb ="org.Mm.eg.db", 
                                pvalueCutoff=0.01,
                                keyType ="SYMBOL")

write.csv(T_NK_MVC_up_GO@compareClusterResult,"T_NK/T_NK_MVC_up_GO_Result.csv")

T_NK_MVC_up_GO_Result <- read.csv("T_NK/T_NK_MVC_up_GO.csv",row.names = 1,check.names = F)
T_NK_MVC_up_GO2 <- T_NK_MVC_up_GO 
T_NK_MVC_up_GO2@compareClusterResult <- T_NK_MVC_up_GO_Result

# T_NK_MVC_up_GO_P <- dotplot(T_NK_MVC_up_GO2) + 
#   scale_color_continuous(low='#E24D36',high='black') + 
#   ggtitle("Up-regulated pathways (Model vs.Control)") + xlab(NULL) + ylab("T/NK cell")
# T_NK_MVC_up_GO_P

T_NK_MVC_up_GO_Result <- read.csv("T_NK/T_NK_MVC_up_GO.csv",row.names = 1,check.names = F)

T_NK_MVC_up_GO_P

T_NK_MVC_up_GO_Result$Cluster <- factor(T_NK_MVC_up_GO_Result$Cluster,
                                         levels = c("T_NK_PBMC","T_NK_Liver","T_NK_Spleen"))
T_NK_MVC_up_GO_Result$Description <- factor(T_NK_MVC_up_GO_Result$Description,
                                             levels = rev(unique(T_NK_MVC_up_GO_Result$Description)))
T_NK_MVC_up_GO_P <- ggplot(T_NK_MVC_up_GO_Result,aes(x=Cluster,y=Description)) +
  geom_point(aes(size=Count,color=p.adjust)) +
  theme_bw() + scale_color_continuous(low='#E24D36',high='black') + 
  ggtitle("Up-regulated pathways (Model vs.Control)") +
  scale_size_continuous(name="Count",range = c(1,10))+
  xlab(NULL) + ylab("Mono/T_NK") +
  theme(axis.text.x = element_text(colour = "black",size = 12, vjust =1 ),
        axis.text.y = element_text(colour = "black",size = 12, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),color = "black",size = 12),
        axis.title.y = element_text(angle=90)) +
  scale_y_discrete(labels=function(y) str_wrap(y,width=25)) 
T_NK_MVC_up_GO_P

Macro_PBMC_MVC <- read.csv("Macro/Macro_PBMC_MVC.csv",row.names = 1,check.names = F)
Macro_PBMC_MVC_up <- row.names(Macro_PBMC_MVC[Macro_PBMC_MVC$type == "Up",])
# Macro_PBMC_MVC_down <- row.names(Macro_PBMC_MVC[Macro_PBMC_MVC$type == "Down",])

Macro_Liver_MVC <- read.csv("Macro/Macro_Liver_MVC.csv",row.names = 1,check.names = F)
Macro_Liver_MVC_up <- row.names(Macro_Liver_MVC[Macro_Liver_MVC$type == "Up",])
# Macro_Liver_MVC_down <- row.names(Macro_Liver_MVC[Macro_Liver_MVC$type == "Down",])

Macro_Spleen_MVC <- read.csv("Macro/Macro_Spleen_MVC.csv",row.names = 1,check.names = F)
Macro_Spleen_MVC_up <- row.names(Macro_Spleen_MVC[Macro_Spleen_MVC$type == "Up",])
# Macro_Spleen_MVC_down <- row.names(Macro_Spleen_MVC[Macro_Spleen_MVC$type == "Down",])

ALL_DEG_up <- list(
  Macro_PBMC = Macro_PBMC_MVC_up,
  Macro_Liver = Macro_Liver_MVC_up,
  Macro_Spleen = Macro_Spleen_MVC_up)

Macro_MVC_up_GO <- compareCluster(ALL_DEG_up,
                                  fun="enrichGO",
                                  OrgDb ="org.Mm.eg.db", 
                                  pvalueCutoff=0.01,
                                  keyType ="SYMBOL")
dotplot(Macro_MVC_up_GO) + 
  scale_color_continuous(low='#E24D36',high='black') + 
  ggtitle("Up-regulated pathways (Model vs.Control)") + xlab(NULL)

write.csv(Macro_MVC_up_GO@compareClusterResult,"Macro/Macro_MVC_up_GO_Result.csv")

Macro_MVC_up_GO_Result <- read.csv("Macro/Macro_MVC_up_GO.csv",row.names = 1,check.names = F)

Macro_MVC_up_GO_Result$Cluster <- factor(Macro_MVC_up_GO_Result$Cluster,
                                         levels = c("Macro_PBMC","Macro_Liver","Macro_Spleen"))
Macro_MVC_up_GO_Result$Description <- factor(Macro_MVC_up_GO_Result$Description,
                                             levels = rev(c("immune receptor activity",
                                                            "immunoglobulin binding",
                                                            "ubiquitin protein ligase binding",
                                                            "MHC class I protein binding",
                                                            "extracellular matrix binding",
                                                            "integrin binding",
                                                            "complement receptor activity",
                                                            "Toll-like receptor binding",
                                                            "STAT family protein binding")))
Macro_MVC_up_GO_P <- ggplot(Macro_MVC_up_GO_Result,aes(x=Cluster,y=Description)) +
  geom_point(aes(size=Count,color=p.adjust)) +
  theme_bw() + scale_color_continuous(low='#E24D36',high='black') + 
  # ggtitle("Up-regulated pathways (Model vs.Control)") + 
  scale_size_continuous(name="Count",range = c(1,10))+
  xlab(NULL) + ylab("Mono/Macro")  +
  theme(axis.text.x = element_text(colour = "black",size = 12, vjust =1 ),
        axis.text.y = element_text(colour = "black",size = 12, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),color = "black",size = 12),
        axis.title.y = element_text(angle=90)) +
  scale_y_discrete(labels=function(y) str_wrap(y,width=25))
Macro_MVC_up_GO_P

Neutro_PBMC_MVC <- read.csv("Neutro/Neutro_PBMC_MVC.csv",row.names = 1,check.names = F)
Neutro_PBMC_MVC_up <- row.names(Neutro_PBMC_MVC[Neutro_PBMC_MVC$type == "Up",])
# Neutro_PBMC_MVC_down <- row.names(Neutro_PBMC_MVC[Neutro_PBMC_MVC$type == "Down",])

Neutro_Liver_MVC <- read.csv("Neutro/Neutro_Liver_MVC.csv",row.names = 1,check.names = F)
Neutro_Liver_MVC_up <- row.names(Neutro_Liver_MVC[Neutro_Liver_MVC$type == "Up",])
# Neutro_Liver_MVC_down <- row.names(Neutro_Liver_MVC[Neutro_Liver_MVC$type == "Down",])

Neutro_Spleen_MVC <- read.csv("Neutro/Neutro_Spleen_MVC.csv",row.names = 1,check.names = F)
Neutro_Spleen_MVC_up <- row.names(Neutro_Spleen_MVC[Neutro_Spleen_MVC$type == "Up",])
# Neutro_Spleen_MVC_down <- row.names(Neutro_Spleen_MVC[Neutro_Spleen_MVC$type == "Down",])

ALL_DEG_up <- list(
  Neutro_PBMC = Neutro_PBMC_MVC_up,
  Neutro_Liver = Neutro_Liver_MVC_up,
  Neutro_Spleen = Neutro_Spleen_MVC_up)

Neutro_MVC_up_GO <- compareCluster(ALL_DEG_up,
                                   fun="enrichGO",
                                   OrgDb ="org.Mm.eg.db", 
                                   pvalueCutoff=0.01,
                                   keyType ="SYMBOL")
dotplot(Neutro_MVC_up_GO) + 
  scale_color_continuous(low='#E24D36',high='black') + 
  ggtitle("Up-regulated pathways (Model vs.Control)") + xlab(NULL)

write.csv(Neutro_MVC_up_GO@compareClusterResult,"Neutro/Neutro_MVC_up_GO_Result.csv")

Neutro_MVC_up_GO_Result <- read.csv("Neutro/Neutro_MVC_up_GO.csv",row.names = 1,check.names = T)
# Neutro_MVC_up_GO2 <- Neutro_MVC_up_GO
# Neutro_MVC_up_GO2@compareClusterResult <- Neutro_MVC_up_GO_Result

# Neutro_MVC_up_GO_P <- dotplot(Neutro_MVC_up_GO2) + 
#   scale_color_continuous(low='#E24D36',high='black') + 
#   ggtitle("Up-regulated pathways (Model vs.Control)") + xlab(NULL)
# 
# 
# as.numeric(unique(Neutro_MVC_up_GO2@compareClusterResult$Cluster))

Neutro_MVC_up_GO_Result$Cluster <- factor(Neutro_MVC_up_GO_Result$Cluster,
                                          levels = c("Neutro_PBMC","Neutro_Liver","Neutro_Spleen"))
Neutro_MVC_up_GO_Result$Description <- factor(Neutro_MVC_up_GO_Result$Description,
                                              levels = rev(c("MHC class I protein binding",
                                                             "CD8 receptor binding",
                                                             "TAP1 binding",
                                                             "TAP2 binding",
                                                             "immunoglobulin binding",
                                                             "antigen binding","cytokine activity",
                                                             "tumor necrosis factor receptor binding")))
Neutro_MVC_up_GO_P <- ggplot(Neutro_MVC_up_GO_Result,aes(x=Cluster,y=Description)) +
  geom_point(aes(size=Count,color=p.adjust)) +
  theme_bw() + scale_color_continuous(low='#E24D36',high='black') + 
  scale_size_continuous(name="Count",range = c(1,10))+
  # ggtitle("Up-regulated pathways (Model vs.Control)") + 
  xlab(NULL) + ylab("Neutro")  +
  theme(axis.text.x = element_text(colour = "black",size = 12, vjust =1 ),
        axis.text.y = element_text(colour = "black",size = 12, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),color = "black",size = 12),
        axis.title.y = element_text(angle=90)) +
  scale_y_discrete(labels=function(y) str_wrap(y,width=25))
Neutro_MVC_up_GO_P

T_NK_MVC_up_GO_P / Macro_MVC_up_GO_P / Neutro_MVC_up_GO_P
