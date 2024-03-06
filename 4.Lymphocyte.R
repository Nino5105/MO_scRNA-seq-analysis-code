# 1.1 load packages
library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggsci)
library(RColorBrewer)
library(ggrepel)
library(patchwork)

# 1.2 load  dataset
data <- readRDS("ALL/Lym.RDS")
table(data$subType0523)

DimPlot(data,group.by = "cell_type_filter") + scale_color_d3() + 
  xlab("") + ylab("") + ggtitle("") # 6*5

data$cell_type_filter <- factor(data$cell_type_filter)
data@active.ident <- data$cell_type_filter
table(data$cell_type_filter)

# T_NK subset
T_NK <- subset(data,subset = subType0523 %in% c("Cd4 T_C1_Lef1",
                                                "Cd4 T_C2_Foxp3",
                                                "Cd4 T_C3_Cxcr5",
                                                "Cd4 T_C4_Ifitm1",
                                                "Cd4 T_C5_Cxcr3",
                                                "Cd8 T_C1_Lef1",
                                                "Cd8 T_C2_Gzmk",
                                                "Cd8 T_C3_Mki67",
                                                "Cd8 T_C4_Il7r",
                                                "Cd8 T_C5_Irf8",
                                                "NK cell"))
table(T_NK$cluster)

DimPlot(T_NK,group.by = "subType0523",pt.size = 0.5,label = T,label.size = 5) + 
  scale_color_manual(values = colors) + 
  xlab("") + ylab("") + ggtitle("") + NoLegend()
saveRDS(T_NK,"T_NK.rds")

# B_cell subset

B_cell <- subset(data,subset = subType0523 %in% c("B cell_C1_Zfas1",
                                                  "B cell_C2_Fcer2a",
                                                  "B cell_C3_Iglc1",
                                                  "B cell_C4_Alas2",
                                                  "B cell_C5_Cxcr4",
                                                  "B cell_C6_Mki67"))

table(B_cell$subType0523)
B_cell@active.ident <- B_cell$subType0523
DimPlot(B_cell,group.by = "subType0523",pt.size = 0.5,label = T,label.size = 5) + 
  scale_color_manual(values = colors) + 
  xlab("") + ylab("") + ggtitle("") + NoLegend()
saveRDS(B_cell,"B_cell.rds")

# Plasma subset
Plasma <- subset(data,subset = subType0523 %in% c("Plasma_C1_Ighg2c",
                                                  "Plasma_C2_Ighm",
                                                  "Plasma_C3_Ighg2b",
                                                  "Plasma_C4_Ighg3",
                                                  "Plasma_C5_Ighg1",
                                                  "Plasma_C6_Igha"))
table(Plasma$subType0523)
Plasma@active.ident <- Plasma$subType0523
DimPlot(Plasma,group.by = "subType0523",pt.size = 0.5,label = T,label.size = 5) + 
  scale_color_manual(values = colors) + 
  xlab("") + ylab("") + ggtitle("") + NoLegend()
saveRDS(Plasma,"Plasma.rds")


colors = c("Cd4 T_C1_Lef1" = "#7876B1FF",
           "Cd4 T_C2_Foxp3"= "#E18727FF",
           "Cd4 T_C3_Cxcr5"= "#6F99ADFF",
           "Cd4 T_C4_Ifitm1"="#20854EFF",
           "Cd4 T_C5_Cxcr3"="#BC3C29FF",
           "Cd8 T_C1_Lef1"="#0072B5FF",
           "Cd8 T_C2_Gzmk"="#71D0F5FF",
           "Cd8 T_C3_Mki67"="#46732EFF",
           "Cd8 T_C4_Il7r"="#075149FF",
           "Cd8 T_C5_Irf8"="#7F7F7FFF",
           "NK cell"="#E377C2FF",
           "B cell_C1_Zfas1"="#3C5488FF",
           "B cell_C2_Fcer2a"="#00A087FF",
           "B cell_C3_Iglc1"="#4DBBD5FF",
           "B cell_C4_Alas2"="#E64B35FF",
           "B cell_C5_Cxcr4"="#FFDC91FF",
           "B cell_C6_Mki67"="#EE4C97FF",
           "Plasma_C1_Ighg2c"="#8C564BFF",
           "Plasma_C2_Ighm"="#9467BDFF",
           "Plasma_C3_Ighg2b"="#D62728FF",
           "Plasma_C4_Ighg3"="#2CA02CFF",
           "Plasma_C5_Ighg1"="#FF7F0EFF",
           "Plasma_C6_Igha"="#1F77B4FF")


# 1.3 DEG analysis 

data <- readRDS("Lym-subType-Type-DEGs.RDS")
head(data)

data$subType <- factor(data$subType,levels = c("Cd4 T_C1_Lef1","Cd4 T_C2_Foxp3","Cd4 T_C3_Cxcr5","Cd4 T_C4_Ifitm1",
                                                               "Cd4 T_C5_Cxcr3","Cd8 T_C1_Lef1","Cd8 T_C2_Gzmk","Cd8 T_C3_Mki67",
                                                               "Cd8 T_C4_Il7r","Cd8 T_C5_Irf8","NK cell","B cell_C1_Zfas1",
                                                               "B cell_C2_Fcer2a","B cell_C3_Iglc1","B cell_C4_Alas2","B cell_C5_Cxcr4",
                                                               "B cell_C6_Mki67","Plasma_C1_Ighg2c",
                                                               "Plasma_C2_Ighm","Plasma_C3_Ighg2b","Plasma_C4_Ighg3","Plasma_C5_Ighg1","Plasma_C6_Igha"))


data_up <- data[data$type == "Up",]
table(data_up$subType,data_up$Com)

data_down <- data[data$type == "Down",]
table(data_down$subType,data_down$Com)

length(intersect(subset(data, subType=="Cd4 T_C1_Lef1" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Cd4 T_C1_Lef1" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Cd4 T_C2_Foxp3" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Cd4 T_C2_Foxp3" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Cd4 T_C3_Cxcr5" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Cd4 T_C3_Cxcr5" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Cd4 T_C4_Ifitm1" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Cd4 T_C4_Ifitm1" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Cd4 T_C5_Cxcr3" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Cd4 T_C5_Cxcr3" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Cd8 T_C1_Lef1" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Cd8 T_C1_Lef1" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Cd8 T_C2_Gzmk" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Cd8 T_C2_Gzmk" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Cd8 T_C3_Mki67" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Cd8 T_C3_Mki67" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Cd8 T_C4_Il7r" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Cd8 T_C4_Il7r" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Cd8 T_C5_Irf8" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Cd8 T_C5_Irf8" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="NK cell" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="NK cell" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="B cell_C1_Zfas1" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="B cell_C1_Zfas1" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="B cell_C2_Fcer2a" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="B cell_C2_Fcer2a" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="B cell_C3_Iglc1" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="B cell_C3_Iglc1" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="B cell_C4_Alas2" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="B cell_C4_Alas2" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="B cell_C5_Cxcr4" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="B cell_C5_Cxcr4" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="B cell_C6_Mki67" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="B cell_C6_Mki67" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Plasma_C1_Ighg2c" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Plasma_C1_Ighg2c" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Plasma_C2_Ighm" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Plasma_C2_Ighm" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Plasma_C3_Ighg2b" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Plasma_C3_Ighg2b" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Plasma_C4_Ighg3" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Plasma_C4_Ighg3" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Plasma_C5_Ighg1" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Plasma_C5_Ighg1" & Com=="ModvsART" & type=="Up")[,'gene']))
length(intersect(subset(data, subType=="Plasma_C6_Igha" & Com=="ModvsCon" & type=="Up")[,'gene'],subset(data, subType=="Plasma_C6_Igha" & Com=="ModvsART" & type=="Up")[,'gene']))

length(intersect(subset(data, subType=="Cd4 T_C1_Lef1" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Cd4 T_C1_Lef1" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Cd4 T_C2_Foxp3" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Cd4 T_C2_Foxp3" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Cd4 T_C3_Cxcr5" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Cd4 T_C3_Cxcr5" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Cd4 T_C4_Ifitm1" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Cd4 T_C4_Ifitm1" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Cd4 T_C5_Cxcr3" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Cd4 T_C5_Cxcr3" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Cd8 T_C1_Lef1" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Cd8 T_C1_Lef1" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Cd8 T_C2_Gzmk" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Cd8 T_C2_Gzmk" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Cd8 T_C3_Mki67" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Cd8 T_C3_Mki67" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Cd8 T_C4_Il7r" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Cd8 T_C4_Il7r" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Cd8 T_C5_Irf8" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Cd8 T_C5_Irf8" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="NK cell" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="NK cell" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="B cell_C1_Zfas1" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="B cell_C1_Zfas1" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="B cell_C2_Fcer2a" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="B cell_C2_Fcer2a" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="B cell_C3_Iglc1" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="B cell_C3_Iglc1" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="B cell_C4_Alas2" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="B cell_C4_Alas2" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="B cell_C5_Cxcr4" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="B cell_C5_Cxcr4" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="B cell_C6_Mki67" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="B cell_C6_Mki67" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Plasma_C1_Ighg2c" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Plasma_C1_Ighg2c" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Plasma_C2_Ighm" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Plasma_C2_Ighm" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Plasma_C3_Ighg2b" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Plasma_C3_Ighg2b" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Plasma_C4_Ighg3" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Plasma_C4_Ighg3" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Plasma_C5_Ighg1" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Plasma_C5_Ighg1" & Com=="ModvsART" & type=="Down")[,'gene']))
length(intersect(subset(data, subType=="Plasma_C6_Igha" & Com=="ModvsCon" & type=="Down")[,'gene'],subset(data, subType=="Plasma_C6_Igha" & Com=="ModvsART" & type=="Down")[,'gene']))


# 3.Bubble plots

DEGs <- read.csv("ALL_DEGs.csv")
head(DEGs)
DEGs$Cell_type <- factor(DEGs$Cell_type,levels = c("Cd4 T_C1_Lef1","Cd4 T_C2_Foxp3","Cd4 T_C3_Cxcr5","Cd4 T_C4_Ifitm1",
                                               "Cd4 T_C5_Cxcr3","Cd8 T_C1_Lef1","Cd8 T_C2_Gzmk","Cd8 T_C3_Mki67",
                                               "Cd8 T_C4_Il7r","Cd8 T_C5_Irf8","NK cell","B cell_C1_Zfas1",
                                               "B cell_C2_Fcer2a","B cell_C3_Iglc1","B cell_C4_Alas2","B cell_C5_Cxcr4",
                                               "B cell_C6_Mki67","Plasma_C1_Ighg2c",
                                               "Plasma_C2_Ighm","Plasma_C3_Ighg2b","Plasma_C4_Ighg3","Plasma_C5_Ighg1","Plasma_C6_Igha"))


p1 <- ggplot(DEGs,aes(x=MVC_up.regulated,y=MVA_up.regulated)) +
  geom_point(aes(color=Cell_type,size = overlaping.up.regulated.DEGs),alpha=1) +
  xlab("Model vs Control") + ylab("Model vs ART") +
  theme_bw() + theme(legend.position = c("right")) +
  geom_text_repel(data = DEGs,aes(label = Cell_type),size = 5,
                  segment.color = "black", show.legend = FALSE ) + 
  scale_color_manual(values = c("#7876B1FF","#E18727FF","#6F99ADFF","#20854EFF","#BC3C29FF",
                                "#0072B5FF","#71D0F5FF","#46732EFF","#075149FF","#7F7F7FFF",
                                "#E377C2FF","#3C5488FF","#00A087FF","#4DBBD5FF","#E64B35FF",
                                "#FFDC91FF","#EE4C97FF","#8C564BFF","#9467BDFF","#D62728FF",
                                "#2CA02CFF","#FF7F0EFF","#1F77B4FF"))

p2 <- ggplot(DEGs,aes(x=MVC_down.regulated,y=MVA_down.regulated)) +
  geom_point(aes(color=Cell_type,size = overlaping.down.regulated.DEGs),alpha=1) +
  xlab("Model vs Control") + ylab("Model vs ART") +
  theme_bw()+theme(legend.position = c("none")) +
  geom_text_repel(data = DEGs,aes(label = Cell_type),size = 5,
                  segment.color = "black", show.legend = FALSE ) + 
  scale_color_manual(values = c("#7876B1FF","#E18727FF","#6F99ADFF","#20854EFF","#BC3C29FF",
                                "#0072B5FF","#71D0F5FF","#46732EFF","#075149FF","#7F7F7FFF",
                                "#E377C2FF","#3C5488FF","#00A087FF","#4DBBD5FF","#E64B35FF",
                                "#FFDC91FF","#EE4C97FF","#8C564BFF","#9467BDFF","#D62728FF",
                                "#2CA02CFF","#FF7F0EFF","#1F77B4FF"))
p1 / p2

# 1.4 GO enrichment analysis

# B cell
B_GO <- read.csv("DEGs_Venn/B cell DEGs_GO.csv")
B_GO$log10Pvalue <- -log10(B_GO$pvalue)
head(B_GO)
B_GO$orgin <- factor(B_GO$orgin,levels = rev(c("PBMC","Liver","Spleen")))

ggbarplot(B_GO, x="Description", y="log10Pvalue", fill = "orgin", color = "white",
          palette =  c("PBMC" = '#53A85F',
                       "Liver" = '#F1BB72',
                       "Spleen" = '#E95C59'),
          sort.val = "asc", 
          sort.by.groups=TRUE, 
          x.text.angle=0, 
          xlab = NULL,ylab = "-log10Pvalue") + coord_flip() # 5x10

# T_NK
T_NK_GO <- read.csv("DEGs_Venn/T NK cell DEGs_GO.csv")
T_NK_GO$log10Pvalue <- -log10(T_NK_GO$pvalue)
head(T_NK_GO)
T_NK_GO$orgin <- factor(T_NK_GO$orgin,levels = rev(c("PBMC","Liver","Spleen")))

ggbarplot(T_NK_GO, x="Description", y="log10Pvalue", fill = "orgin", color = "white",
          palette =  c("PBMC" = '#53A85F',
                       "Liver" = '#F1BB72',
                       "Spleen" = '#E95C59'),
          sort.val = "asc", 
          sort.by.groups=TRUE, 
          x.text.angle=0, 
          xlab = NULL,ylab = "-log10Pvalue") + coord_flip() # 8*5

# Plasma
Plasma_GO <- read.csv("DEGs_Venn/Plasma DEGs_GO.csv")
Plasma_GO$log10Pvalue <- -log10(Plasma_GO$pvalue)
head(Plasma_GO)
Plasma_GO$orgin <- factor(Plasma_GO$Orgin,levels = rev(c("Liver","Spleen")))

ggbarplot(Plasma_GO, x="Description", y="log10Pvalue", fill = "orgin", color = "white",
          palette =  c("PBMC" = '#53A85F',
                       "Liver" = '#F1BB72',
                       "Spleen" = '#E95C59'),
          sort.val = "asc", 
          sort.by.groups=TRUE, 
          x.text.angle=0, 
          xlab = NULL,ylab = "-log10Pvalue") + coord_flip() # 8*5