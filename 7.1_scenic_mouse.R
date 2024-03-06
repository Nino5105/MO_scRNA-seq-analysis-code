library(tidyverse)
library(patchwork)
library(SCENIC)
library(foreach)
library(Seurat)
library(pheatmap)
library(argparser)

argv<-arg_parser("")
argv<-add_argument(argv,"--rds",help="expr_path")
#argv<-add_argument(argv,"--meta",help="meta_csv_path")
argv<-add_argument(argv,"--outdir",help="the output dir")
argv<-parse_args(argv)

##transfer parameters

input <- argv$rds
exprMat<-readRDS(input)

print(dim(exprMat))
print(exprMat[1:4,1:4])
#print(head(scRNA@meta.data))
#print(dim(scRNA))

outdir<-argv$outdir

dir.create(paste0(outdir,'/',"SCENIC"))
dir.create(paste0(outdir,'/',"SCENIC/int"))
setwd(paste0(outdir,'/',"SCENIC"))

# meta<-argv$meta

# cellInfo <- read.csv(meta)
# colnames(cellInfo)[which(colnames(cellInfo)=="stim")] <- "sample"
# #colnames(cellInfo)[which(colnames(cellInfo)=="class")] <- "class"
# colnames(cellInfo)[which(colnames(cellInfo)=="cell")] <- "celltype"
# #cellInfo <- cellInfo[,c("sample","class","celltype")]
# print(head(cellInfo))
# saveRDS(cellInfo, file="int/cellInfo.Rds")
#exprMat <- as.matrix(scRNA@assays$RNA@counts)

mydbDIR <- "/szrmyy/wangjgLab/scRNA/wangc/brain_Rat/ref/cisTarget_databases"
mydbs <- c("mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
           "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")

scenicOptions <- initializeScenic(org="mgi",
                                  nCores=8,
                                  dbDir=mydbDIR,
                                  dbs = mydbs,
                                  datasetTitle = "AA")

print("creat scenicOptions!")

genesKept <- geneFiltering(exprMat, scenicOptions,
              minCountsPerGene = 3 * 0.01 * ncol(exprMat),
              minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)
save(exprMat_filtered,scenicOptions,file = "input_GENIE3_data.Rdata")

print("runGenie3 done!")
