#Pimchanok Phankeaw
#whole transcriptomic sequencing from tissue samples containing various cell types
#Compare DEGs between two groups

#load library
library(Seurat)
library(ggplot2)
library(stringr)
library(openxlsx)

#set working directory
setwd('...')
#load gene counts
data <- read.table('....txt',header=TRUE,sep='\t')

#create seurat object
mydata <- CreateSeuratObject(counts = data, min.cells = 0, min.features = 20, project = '...')
mydata <- NormalizeData(mydata,normalization.method = 'CLR',scale.factor = 10000,margin = 1, block.size = NULL,verbose = TRUE)
mydata <- FindVariableFeatures(mydata, selection.method="mean.var.plot")
mydata <- ScaleData(mydata,features=VariableFeatures(mydata))

#read meta data
meta_data <- read.table('....txt,)

#add sample info to seurat object
mydata$... <- as.character(meta_data[1,])
mydata$... <- as.character(meta_data[2,])
mydata$... <- as.character(meta_data[3,])

#GET FRACTION/COMPOSITION VALUES 
fractiondataset <- as.data.frame(table(mydata@meta.data$..., mydata@meta.data$..., mydata@meta.data$...))

#start code here

#ESTABLISH GROUPING
DefaultAssay(mydata) <- "RNA"
#GROUP1 cell1
Idents(mydata) <- "cell"
cell1 <- subset(mydata, idents = "cell1")
#GROUP2 fibroblasts
Idents(mydata) <- "cell"
cell2 <- subset(mydata, idents = "cell2")

#next set ident each cell in each group
#GROUP3 cell1 in group1
Idents(cell1) <- "condition1"
cell1group1 <- subset(cell1, idents = "group1")
#GROUP4 cell1 in group2
Idents(cell1) <- "condition1"
cell1group2 <- subset(cell1, idents = "group2")

#Differential expression analysis
#comparison 1 group1 vs group2 
Idents(cell1) <- "condition1"
group1vsgroup2 <- FindMarkers(cell1, ident.1 = 'group1', ident.2 = 'group2', test.use = "DESeq2", only.pos = FALSE,min.pct = -Inf, logfc.threshold = -Inf, min.cells.feature = 1, min.cells.group = 1,assay= "RNA")

write.xlsx(group1vsgroup2,'your directory/filesname.xlsx,col.names = TRUE, row.names = TRUE')
