rm(list=ls())
library(Seurat)
library(enrichplot)
#library(tidyverse)
#library(ggstatsplot)
setwd("/data01/Rawdata/data/wangchen/multiorgan/gzmk")

rds<-readRDS("gzmb.rds")
head(rds@meta.data)
unique(rds$type)
Idents(rds)<-rds$type
diffgenes<-FindMarkers(rds,ident.1 = "ART", ident.2 = "Model")
head(diffgenes)
diffgenes<-subset(diffgenes,p_val_adj<=0.05)
diffgenes<-subset(diffgenes,avg_log2FC >= 0.25)
head(diffgenes)
VlnPlot(rds,features = "H2-K1",group.by ="type")

library(clusterProfiler)
library(org.Mm.eg.db)

trans<-bitr(rownames(diffgenes),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db")

ekk <- enrichKEGG(gene = trans$ENTREZID,
                  keyType = 'kegg',
                  organism = 'mmu',
                  pvalueCutoff = 0.05,
                  pAdjustMethod  = "BH"
                  #qvalueCutoff  = 0.05
)

# 一句代码转化，kegg为富集结果
kegg_SYMBOL <- DOSE::setReadable(ekk,
                                 OrgDb = "org.Mm.eg.db",
                                 keyType = "ENTREZID")

head(kegg_SYMBOL@result)
write.csv(kegg_SYMBOL@result,"gzmk_artvsmodel_kegg.csv")
dotplot(kegg_SYMBOL,showCategory =20)
#
organism = 'mmu'    #  人类'hsa' 小鼠'mmu'   
OrgDb = 'org.Mm.eg.db'#人类"org.Hs.eg.db" 小鼠"org.Mm.eg.db"

diffgenes<-FindMarkers(rds,ident.1 = "ART", ident.2 = "Model")
diffgenes<-subset(diffgenes,p_val_adj<=0.05)
trans<-bitr(rownames(diffgenes),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db")

diffgenes$SYMBOL<-rownames(diffgenes)

need_DEG <- merge(diffgenes, trans, by='SYMBOL')  #按照SYMBOL合并注释信息
head(need_DEG)
geneList <- need_DEG$avg_log2FC
names(geneList) <- need_DEG$ENTREZID
geneList <- sort(geneList, decreasing = T)
KEGG_kk_entrez <- gseKEGG(geneList     = geneList,
                          organism     = organism, #人hsa 鼠mmu
                          pvalueCutoff = 0.25)  #实际为padj阈值可调整 
head(KEGG_kk_entrez)
KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez, 
                             OrgDb=OrgDb,
                             keyType='ENTREZID')#转化id   
##选取富集结果
kk_gse <- KEGG_kk
kk_gse_entrez <- KEGG_kk_entrez

###条件筛选 
#一般认为|NES|>1，NOM pvalue<0.05，FDR（padj）<0.25的通路是显著富集的
kk_gse_cut <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$p.adjust<0.25 & abs(kk_gse$NES)>1]
kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0,]
kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0,]

#选择展现NES前几个通路 
#down_gsea <- kk_gse_cut_down[tail(order(kk_gse_cut_down$NES,decreasing = T),10),]
#up_gsea <- kk_gse_cut_up[head(order(kk_gse_cut_up$NES,decreasing = T),10),]
#diff_gsea <- kk_gse_cut[head(order(abs(kk_gse_cut$NES),decreasing = T),10),]
pathways<-c("mmu04630","mmu04668")
down_gsea <- kk_gse_cut_down[kk_gse_cut_down$ID %in% pathways,]
#

#devtools::install_github( "junjunlab/GseaVis") 

library(stringr)
library(clusterProfiler)
library(GseaVis)

genes<-c("Crlf2","Stat3","Socs3","Ifngr1","Raf1","Cdkn1a","Pim1")

pdf("gsea_down.pdf")
gseaNb(object = kk_gse,
       geneSetID = down_gsea$Description,
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8),
       addGene = genes)
dev.off()

gseaNb(object = kk_gse,
       geneSetID = "JAK-STAT signaling pathway - Mus musculus (house mouse)",
       subPlot = 2,
       addGene = genes)

saveRDS(kk_gse,"kk_gse_gzmb.rds")
gseap2 <- gseaplot2(kk_gse,
                    down_gsea$ID,#富集的ID编号
                    title = "DN_GSEA_all",#标题
                    color = "blue",#GSEA线条颜色
                    base_size = 20,#基础字体大小
                    rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                    subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                    ES_geom = "line",#enrichment score用线还是用点"dot"
                    pvalue_table = T) #显示pvalue等信息

ggsave(gseap2, filename = "GSEA_ARTvsModel.pdf",width =10,height =7)
write.csv(kk_gse,"gzmk_kegg_gse.csv")
gene<-c("Stat3","Il10","Crlf2","Socs3","Ifngr1","Il2ra","Il2rb","Raf1","Cdkn1a","Pim1")
VlnPlot(rds,features = gene,group.by ="type")


#### 经典的GSEA图 
up_gsea$Description
i=2
gseap1 <- gseaplot2(kk_gse,
                    up_gsea$ID[i],#富集的ID编号
                    title = up_gsea$Description[i],#标题
                    color = "red", #GSEA线条颜色
                    base_size = 20,#基础字体大小
                    rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                    subplots = 1:3,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                    ES_geom = "line", #enrichment score用线还是用点"dot"
                    pvalue_table = T) #显示pvalue等信息
ggsave(gseap1, filename = 'GSEA_up_1.pdf', width =10, height =8)

#### 合并 GSEA通路 
gseap2 <- gseaplot2(kk_gse,
                    up_gsea$ID,#富集的ID编号
                    title = "UP_GSEA_all",#标题
                    color = "red",#GSEA线条颜色
                    base_size = 20,#基础字体大小
                    rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                    subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                    ES_geom = "line",#enrichment score用线还是用点"dot"
                    pvalue_table = T) #显示pvalue等信息
ggsave(gseap2, filename = "GSEA_up_all.pdf",width =12,height =12)

###########
head(rds@meta.data)
str(rds)
?AverageExpression
Idents(rds)<-rds$orgin
unique(rds$orgin)
sce<-subset(rds,idents="Liver")
av <-AverageExpression(sce,assays = "RNA", group.by ="sample")
av<-av$RNA
av_gene<-av["Ccl3",]
av_gene_sc<-data.frame(RNA_seq=av["Ccl3",],sample=colnames(av))
head(av_gene)
unique(av_gene$sample)
av_gene[av_gene$sample %in% c("LA1","LA2","LA3")]<-

  dt1 <- dt[,c(1,3,4)]
dt2 <- dt[,1:2]
#查看数据；
dt1
dt2
#载入相关的R包；
library(ggplot2)
library(ggpubr)
library(cowplot)

#使用数据1绘制柱状图；
p1 <- ggplot(dt1, aes(x = sample, y = Q_PCR,ymin=Q_PCR-SE,ymax=Q_PCR+SE))+
  geom_bar(stat = "identity", fill = "#FF99CC99",
           color="#FF99CC",width = 0.6)+
  scale_x_discrete(name = NULL)+
  geom_errorbar(width = 0.2)+
  theme_classic()
p1

p2 <- p1 + scale_y_continuous(expand = expansion(mult = c(0, 0)),
                              limits = c(0,70),
                              breaks = c(0,10,20,30,40,50,60,70),
                              label = c("0","10","20","30","40","50","60","70"),
                              position = "left",name = "Q-PCR relative expression")
p2

p3 <- ggplot(data=dt2, aes(x = sample, y = RNA_seq)) +
  geom_line(aes(group = 1),color="#c77cff",
            linewidth=1,show.legend=T) +
  geom_point(color="#c77cff",shape=21,fill="white",
             size=3,show.legend=F)
p3

p4 <- p3 + scale_y_continuous(expand = expansion(mult = c(0, 0)),
                              limits=c(0, 700),
                              breaks = c(0,100,200,300,400,500,600,700),
                              label = c("0","100","200","300","400","500","600","700"),
                              position = "right",name = "RNA-seq expression") +
  scale_x_discrete(name = NULL)+
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks")
p4

two_plots <- align_plots(p2, p4, align="hv", axis="tblr")
#使用draw_plot在图层上添加另一个图层；
p5 <- ggdraw(two_plots[[1]]) + draw_plot(two_plots[[2]])
p5

#自定义图表主题，对图表做精细调整；
top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2
#绘图区域留白；
mytheme<-theme(plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                                units="inches"))
p5+mytheme






