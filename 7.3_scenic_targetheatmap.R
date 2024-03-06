rds<-readRDS("/public/Count/scRNA/multi_organ/Scenic/Tcell/Cd8_C2_Gzmk_exprMat_DEG/SCENIC/int/2.5_regulonTargetsInfo.Rds")
head(rds)
unique(rds$TF)

tf_pick<-rds[rds$TF %in% c("Stat3","Foxp1"),]
tf_pick<-subset(tf_pick,highConfAnnot=="TRUE")

write.csv(tf_pick,"/public/Count/scRNA/multi_organ/6.gzmb/tf_pick.csv")

gzmb<-readRDS("/public/Count/scRNA/multi_organ/6.gzmb/gzmb.rds")

ery.avg = AverageExpression(gzmb, assays='SCT', group.by ="GroupAll", slot = "data")
ery.avg = as.data.frame(ery.avg$SCT)

tf_pick_1<-tf_pick[tf_pick$TF %in% c("Stat3"),]
heatmap_data<-ery.avg[tf_pick_1$gene,]

pheatmap(heatmap_data,scale="row",cluster_col=F,main="Stat3")


tf_pick_1<-tf_pick[tf_pick$TF %in% c("Foxp1"),]
heatmap_data<-ery.avg[tf_pick_1$gene,]

pheatmap(heatmap_data,scale="row",cluster_col=F,main="Foxp1")

