multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
cols = c(
  "#F0A3FF", #Amethyst
  "#0075DC", #Blue
  "#993F00", #Caramel
  "#4C005C", #Damson
  "#191919", #Ebony
  "#005C31", #Forest
  "#2BCE48", #Green
  "#FFCC99", #Honeydew
  "#808080", #Iron
  "#94FFB5", #Jade
  "#8F7C00", #Khaki
  "#9DCC00", #Lime
  "#C20088", #Mallow
  "#003380", #Navy
  "#FFA405", #Orpiment
  "#FFA8BB", #Pink
  "#426600", #Quagmire
  "#FF0010", #Red
  "#5EF1F2", #Sky
  "#00998F", #Turquoise
  "#E0FF66", #Uranium
  "#740AFF", #Violet
  "#990000", #Wine
  "#FFFF80", #Xanthin
  "#FFFF00", #Yellow
  "#FF5005" #Zinnia
)
require(ggplot2)
theme_set(theme_bw())
require(ggrepel)
require(stringr)
require(reshape2)
require(ggtree)
require(ape)
require(phytools)
## load libraries for the heat map
library("RColorBrewer")
library("gplots")
library("vsn")
require("pheatmap")
library("AnnotationDbi")
library("genefilter")
library("topGO")
require("DESeq2")
require("VennDiagram")
require("Hmisc")
require("PopGenome")
require("NOISeq") ###-- Diagnostics plot
require(cqn) ###--- Normalization for GC content and gene length if gene x sample bias

#############--- ANALYSIS RNA SEQ DATA ---##############
draft = '~/Documents/INRA/GEMANEMA/Paper/'
setwd('~/Documents/INRA/GEMANEMA/Paper/Data/RNA/')

library("tximport")
library("readr")
library("tximportData")
library("GenomicFeatures")
library("DESeq2")
library("ggplot2")

#------=====#### Sample information + Gene to transcript association ####=====------
sampInfo0 = read.table('samplelist.txt',header=F,sep='\t')
colnames(sampInfo0) = c('samp')
sampInfo0$sample = paste0(sapply(stringr::str_split(sampInfo0$samp,"_"),function(x) x[2]),row.names(sampInfo0))
files <- file.path("./RNAsalm",paste0(sampInfo0$sample,"_quant"), "quant.sf")
all(file.exists(files))

sampInfo0$sheep = sapply(stringr::str_split(sampInfo0$samp,"_"),function(x) x[2])

extract = read.table(file='../DNA_RNA_extraction.txt',header=T)
sampInfo0$condition = extract$Group[match(sampInfo0$sheep, extract$SheepID)]
sampInfo0$rep = sapply(str_split(sampInfo0$samp,"_"),function(x) x[3])

## Rename 432 into 232 (mislabelled)
#sampInfo0$sheep[sampInfo0$sheep==432]=232

colsRS = c('#e08214','#8073ac',
           '#d73027','#4575b4')
colEnv = c('#bf812d','#35978f',
           '#c51b7d','#7fbc41') 

##--- Correspondance between Transcript to Gene ID
gffFile <- "~/Documents/REF_FILES/haemonchus_contortus.PRJEB506.WBPS13.annotations.gff3.gz"
txdb <- makeTxDbFromGFF(file = gffFile,
                        dataSource="WBSP13",
                        organism="Haemonchus contortus")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME") ## dplyr has the same function

#------=====#### DESeq2 / SALMON analysis - gene-level analysis ####=====------

##-- 2291_quant does not behave like other samples (38% variance): removed from further analyses
sampInfo = sampInfo0[-1,]
files <- file.path("./RNAsalm",paste0(sampInfo$sample,"_quant"), "quant.sf")
all(file.exists(files))

##-- import Salmon files
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)

##-- Create dds object
ddsTxi <- DESeq2::DESeqDataSetFromTximport(txi,
                                   colData = sampInfo,
                                   design = ~ condition)
ddsTxi
# class: DESeqDataSet 
# dim: 19384 23 
# metadata(1): version
# assays(2): counts avgTxLength
# rownames(19384): HCON_00000010 HCON_00000020 ... HCON_00194310 HCON_00665875
# rowData names(0):
# colnames: NULL
# colData names(5): samp sample sheep condition rep

###--- Filter counts at least 5 counts in 3 samples
dds <- estimateSizeFactors(ddsTxi)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5)>= 3
dds <- ddsTxi[idx,]

rld <- rlog(ddsTxi)
data <- plotPCA(rld, intgroup = c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
data$name = paste(sampInfo$condition,sampInfo$sheep,sampInfo$rep,sep="_")

ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) +
  geom_text(aes(label=data$name),vjust=1.5,hjust=0.5)+
  theme_bw() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position=c(0,0),
        title=element_text(size=18)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))


####-- Heatmap of the sample-to-sample distances
sampleDists <- dist(t( assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )

rownames(sampleDistMatrix) <- paste(sampInfo$condition,sampInfo$sheep,sampInfo$rep,sep="_")
colnames(sampleDistMatrix) <- paste(sampInfo$condition,sampInfo$sheep,sampInfo$rep,sep="_")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,cex=1.5)

###--- Estimate DE
# run DEseq
dds.salm <- DESeq(dds)

# get differentially expressed genes
resultsNames(dds.salm)
#[1] "Intercept"        "condition_S_vs_R"  /// R are reference

res05.salm <- results(dds.salm, alpha=.05)
table(res05.salm$padj < .05)
# FALSE  TRUE 
# 12456    20 

summary(res05.salm)
# out of 12476 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 7, 0.056%
# LFC < 0 (down)     : 13, 0.1%
# outliers [1]       : 0, 0%
# low counts [2]     : 1, 0.008%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

resOrdered.salm <- res05.salm[order(res05.salm$padj),]
sig.salm = resOrdered.salm[!is.na(resOrdered.salm$padj) &
                   resOrdered.salm$padj<0.05 &
                   abs(resOrdered.salm$log2FoldChange)>=0,]

sigFC1.salm <- resOrdered.salm[!is.na(resOrdered.salm$padj) &
                       resOrdered.salm$padj<0.05 &
                       abs(resOrdered.salm$log2FoldChange)>=1,]
dim(sigFC1.salm)
#[1] 18  6

sigFC2.salm <- resOrdered.salm[!is.na(resOrdered.salm$padj) &
                       resOrdered.salm$padj<0.05 &
                       abs(resOrdered.salm$log2FoldChange)>=2,]
dim(sigFC2.salm)
# [1] 7  6

# ###--- Any bias here ?
# plotMA(resOrdered.salm, ylim=c(-5,5),main="MAplot S vs R sheep")
# legend('topright','Log2FC > 2',pch=1,col='dodgerblue')
# #topGene <- rownames(sigFC2)
# topGene.salm <- seq(1:nrow(sigFC2.salm))
# with(sigFC2.salm[topGene.salm, ], {
#   points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
#   text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
# })

###--- GOI Salmon
selected.salm = factor(rownames(sigFC1.salm))
selected.salm
# [1] HCON_00007274 HCON_00087240 HCON_00182130 HCON_00132830 HCON_00082860 HCON_00050580 HCON_00016780
# [8] HCON_00125260 HCON_00029050 HCON_00109900 HCON_00053620 HCON_00191700 HCON_00007272 HCON_00090510
# [15] HCON_00006490 HCON_00039260 HCON_00063030 HCON_00090580

topsig.salm = data.frame(t(txi$abundance[which(row.names(txi$abundance) %in% selected.salm),]))

topsig.salm$Group = sampInfo$condition
topsig.salm$samp = sampInfo$sheep
m1.salm = melt(topsig.salm,c('Group','samp'))

ggplot(m1.salm,aes(x=variable,y=value,color=Group)) +
  geom_boxplot() +
  scale_y_log10() +
  coord_flip()
  # facet_wrap(~ variable, nrow=3)

###--- Any change in isoform usage ? -> none
tx2gene[tx2gene$GENEID %in% selected.salm,]
#              TXNAME        GENEID
# 1543  HCON_00029050 HCON_00029050
# 2023  HCON_00191700 HCON_00191700
# 2219  HCON_00192150 HCON_00192150
# 2724  HCON_00016780 HCON_00016780
# 4765  HCON_00053620 HCON_00053620
# 5703  HCON_00039260 HCON_00039260
# 6304  HCON_00050580 HCON_00050580
# 8077  HCON_00082580 HCON_00082580
# 8089  HCON_00082860 HCON_00082860
# 8330  HCON_00087240 HCON_00087240
# 9289  HCON_00074250 HCON_00074250
# 10175 HCON_00090580 HCON_00090580
# 13937 HCON_00125260 HCON_00125260
# 14996 HCON_00144960 HCON_00144960
# 16117 HCON_00132830 HCON_00132830
# 19179 HCON_00190640 HCON_00190640
# 20211 HCON_00182130 HCON_00182130
rm(ddsTxi, dds.salm)

#------=====#### STAR counts analysis // DESeq2 #####=====------

###--- Format geneCount table keeping lane information
countData = read.table(file=paste("./RNAStarMulti/229_1.ReadsPerGene.out.tab",sep=""),sep="\t",header=F)
countData = countData[-c(1:4),c(1,4)]
colnames(countData)=c('geneID','229_1')

sampList = paste0(sampInfo0$sheep,"_",rownames(sampInfo0))
sampList = sampList[-1]

#- Format gene counts
for(i in sampList){
  tmp = read.table(file=paste("./RNAStarMulti/",i,".ReadsPerGene.out.tab",sep=""),sep="\t",header=F)
  tmp = tmp[-c(1:4),c(1,4)]
  colnames(tmp) = c('geneID',i)
  countData = merge(countData,tmp,by='geneID')
  rm(tmp)
}
sheepList = factor(sampInfo$sheep)
rownames(countData) = countData$geneID
countData$geneID = NULL

## Remove 0 counts
countData = countData[which(rowSums(countData)!=0),]

##- Get rid of 229_1 / outlier on PCA
countData = countData[,-1]

###------- QC control overall study---------#####

geneinfo = read.table('./gene_information_WBP13.txt',header=T,sep='\t')
geneinfo = geneinfo[geneinfo$Gene.stable.ID %in% rownames(countData),]
mylength = (abs(geneinfo$Gene.end..bp.-geneinfo$Gene.start..bp.))
names(mylength) = geneinfo$Gene.stable.ID
mygc = geneinfo$X..GC.content
names(mygc) = geneinfo$Gene.stable.ID
mydata <- readData(data = countData, length = mylength, gc = mygc,factors=sampInfo)

##--- Saturation
mysaturation = dat(mydata, k = 0, ndepth = 7, type = "saturation")
par(mfrow=c(2,2))
explo.plot(mysaturation, toplot = 1, samples = 1:6, yleftlim = NULL, yrightlim = NULL)
explo.plot(mysaturation, toplot = 1, samples = 7:12, yleftlim = NULL, yrightlim = NULL)
explo.plot(mysaturation, toplot = 1, samples = 13:18, yleftlim = NULL, yrightlim = NULL)
explo.plot(mysaturation, toplot = 1, samples = 19:24, yleftlim = NULL, yrightlim = NULL)

##--- Sensitivity plot
mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
par(mfrow=c(1,1))
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")
mycountsbio2 = dat(mydata, factor = "condition", type = "countsbio")
explo.plot(mycountsbio2, toplot = 1, samples = NULL, plottype = "barplot")

##--- Length bias by condition
mylengthbias = dat(mydata, factor = NULL, type = "lengthbias")
par(mfrow=c(2,2))
explo.plot(mylengthbias, samples = 1:6, toplot = "global")
explo.plot(mylengthbias, samples = 7:12, toplot = "global")
explo.plot(mylengthbias, samples = 13:18, toplot = "global")
explo.plot(mylengthbias, samples = 18:24, toplot = "global")

###--- Length bias varies across samples / but not by condition

# ##- Length bias by sheep
# mylengthbias = dat(mydata, factor = "sheepID", type = "lengthbias")
# explo.plot(mylengthbias, samples = NULL, toplot = "global")
#
# ##- Length bias by nworm
# mylengthbias = dat(mydata, factor = "nworm", type = "lengthbias")
# explo.plot(mylengthbias, samples = NULL, toplot = "global")
#
# ##- Length bias by condition
# mylengthbias = dat(mydata, factor = "condition", type = "lengthbias")
# explo.plot(mylengthbias, samples = NULL, toplot = "global")

###--- GC content bias
myGCbias = dat(mydata, factor = NULL, type = "GCbias")
par(mfrow=c(2,2))
explo.plot(myGCbias, samples = 1:6, toplot = "global")
explo.plot(myGCbias, samples = 7:12, toplot = "global")
explo.plot(myGCbias, samples = 13:18, toplot = "global")
explo.plot(myGCbias, samples = 18:24, toplot = "global")

##-- GC bias looks rather the same across samples but for Hc10_429_1

#####-- PCA
myPCA = dat(mydata, type = "PCA")
explo.plot(myPCA, factor = "condition")

####========================================================
####----- Perform normalization for GC and gene length
####========================================================
###-- From now on, discard genes with no counts
# 
# uCovar=cbind(mylength,mygc)
# uCovar=data.frame(uCovar)
# uCovar = uCovar[sort(rownames(uCovar)),]
# uCovar = uCovar[which(rownames(uCovar) %in% rownames(countData)),]
# countData = countData[which(rownames(countData) %in% rownames(uCovar)),]
# countData = countData[sort(rownames(countData)),]
# sizeFactors = colSums(countData)
# names(sizeFactors)=colnames(countData)
# 
# stopifnot(all(rownames(countData) == rownames(uCovar)))
# stopifnot(colnames(countData) == names(sizeFactors))
# 
# cqnObject <- cqn(countData, lengths = uCovar$mylength,
#                   x = uCovar$mygc, sizeFactors = sizeFactors,
#                   verbose = TRUE)
# cqnOffset <- cqnObject$glm.offset
# cqnNormFactors <- exp(cqnOffset)
# 
# ###--- DESeq2 ---#####
# 
# #------- create DESeq input matrix
# dds0 <- DESeqDataSetFromMatrix(countData, sampInfo,
#                                formula(~ condition))
# nrow(dds0) # 18207
# 
# ##------- Apply cqn normalization
# normFactors <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))
# normalizationFactors(dds) <- normFactors
# 
# ##--- Transformed values for visualization --> Hc10_229_1 to be discarded
# rld <- rlog(dds0)
# data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
# percentVar <- round(100 * attr(data, "percentVar"))
# ggplot(data, aes(PC1, PC2, color=condition)) +
#   geom_point(size=4) +
#   geom_text(aes(label=data$name),vjust=1.5,hjust=0.5)+
#   theme_bw() +
#   theme(axis.text=element_text(size=16),
#         axis.title=element_text(size=16),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=14),
#         legend.position=c(0,0),
#         title=element_text(size=18)) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance"))
# rm(rld,dds0)

###===================================================================
#####------------- ANALYSIS RUN without 229_1
###===================================================================
# 
# countData = countData[,-1]
# sampInfo = sampInfo[-1,]

uCovar = cbind(mylength,mygc)
uCovar = data.frame(uCovar)
uCovar = uCovar[sort(rownames(uCovar)),]
uCovar = uCovar[which(rownames(uCovar) %in% rownames(countData)),]
sizeFactors = colSums(countData)
names(sizeFactors)=colnames(countData)

stopifnot(all(rownames(countData) == rownames(uCovar)))
stopifnot(colnames(countData) == names(sizeFactors))

cqnObject <- cqn::cqn(countData, lengths = uCovar$mylength,
                 x = uCovar$mygc, sizeFactors = sizeFactors,
                 verbose = TRUE)
cqnOffset <- cqnObject$glm.offset
cqnNormFactors <- exp(cqnOffset)

###--- DESeq2 ---#####

#------- create DESeq input matrix
dds <- DESeqDataSetFromMatrix(countData, sampInfo,
                               formula(~ condition))
nrow(dds) 
#[1] 18238

##------- Apply cqn normalization
normFactors <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))
normalizationFactors(dds) <- normFactors

rld <- rlog(dds)
data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) +
  geom_text(aes(label=data$name),vjust=1.5,hjust=0.5)+
  theme_bw() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position=c(0,0),
        title=element_text(size=18)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

##------ Check what's going on with Hc20_422 sample
par(mfrow=c(2,1))
explo.plot(myGCbias, samples = c(13:23), toplot = "global") ##-- 422_16 behave normally
explo.plot(mylengthbias, samples = c(13:23), toplot = "global")

##------- Analysis without merging replicates 
dds0 <- DESeqDataSetFromMatrix(countData, sampInfo,
                              formula(~ condition))
nrow(dds0) 
#[1] 18238

###--- Filter counts at least 5 counts in 3 samples
dds0 <- estimateSizeFactors(dds0)
idx <- rowSums( counts(dds0, normalized=TRUE) >= 5)>= 3
dds <- dds0[idx,]

dds
# class: DESeqDataSet 
# dim: 13561 23 
# metadata(1): version
# assays(1): counts
# rownames(13561): HCON_00000010 HCON_00000020 ... HCON_00194310 HCON_00665875
# rowData names(0):
# colnames(23): 229_2 229_3 ... 432_23 432_24
# colData names(6): samp sample ... rep sizeFactor

####-- Heatmap of the sample-to-sample distances // no particular clustering associated with our condition
sampleDists <- dist(t( assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(sampInfo$condition,sampInfo$sheep,sep="_")
colnames(sampleDistMatrix) <- paste(sampInfo$condition,sampInfo$sheep,sep="_")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,cex=1.2)

###--- Estimate DE
# run DEseq
dds <- DESeq(dds)

# get differentially expressed genes
resultsNames(dds)
#[1] "Intercept"        "condition_S_vs_R"  /// R are reference

res05 <- results(dds, alpha=.05)
table(res05$padj < .05)
# FALSE  TRUE 
# 13428   198 

summary(res05)
# out of 13561 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 57, 0.42%
# LFC < 0 (down)     : 142, 1%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

resOrdered <- res05[order(res05$padj),]
sig05.star = resOrdered[!is.na(resOrdered$padj) &
                   resOrdered$padj<0.05 &
                   abs(resOrdered$log2FoldChange)>=0,]

sigFC1 <- resOrdered[!is.na(resOrdered$padj) &
                       resOrdered$padj<0.05 &
                       abs(resOrdered$log2FoldChange)>=1,]
dim(sigFC1)
#[1] 38  6
summary(sigFC1)
# adjusted p-value < 0.05
# LFC > 0 (up)       : 9, 24%
# LFC < 0 (down)     : 29, 76%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

sigFC2 <- resOrdered[!is.na(resOrdered$padj) &
                       resOrdered$padj<0.05 &
                       abs(resOrdered$log2FoldChange)>=2,]
dim(sigFC2)
# [1] 7  6

###--- Any bias here ?
DESeq2::plotMA(resOrdered, main="MAplot S vs R sheep")
legend('topright','Log2FC > 2',pch=1,col='dodgerblue')
topGene <- seq(1:nrow(sigFC2))
with(sigFC2[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

###--- Results
selected.star = factor(rownames(sigFC1))
selected.star
# [1] HCON_00132830 HCON_00191400 HCON_00007272 HCON_00182130 HCON_00193080 HCON_00007274 HCON_00192750 HCON_00087240 HCON_00090510
# [10] HCON_00050580 HCON_00121620 HCON_00193430 HCON_00134140 HCON_00029050 HCON_00106640 HCON_00157290 HCON_00090580 HCON_00115500
# [19] HCON_00191530 HCON_00191700 HCON_00035530 HCON_00008680 HCON_00048340 HCON_00017230 HCON_00109330 HCON_00035210 HCON_00098360
# [28] HCON_00074200 HCON_00090630 HCON_00111730 HCON_00008690 HCON_00184620 HCON_00008580 HCON_00008540 HCON_00008660 HCON_00039260
# [37] HCON_00193200 HCON_00036580
topsig.star = data.frame(t(countData[which(rownames(countData) %in% selected.star),]))

topsig.star$Group = sampInfo$condition
topsig.star$samp = rownames(topsig.star)
m1 = melt(topsig.star,c('Group','samp'))

ggplot(m1,aes(x= variable, y = value,col = Group)) +
  geom_boxplot() + coord_flip() + scale_y_log10()

#------=====#### STAR counts analysis // VOOM #####=====------
detach("package:DESeq2", unload = T)
require("limma")
require("edgeR")

## import count data
countData = read.table(file=paste("./RNAStarMulti/229_1.ReadsPerGene.out.tab",sep=""),sep="\t",header=F)
countData = countData[-c(1:4),c(1,4)]
colnames(countData)=c('geneID','229_1')

sampList = paste0(sampInfo0$sheep,"_",rownames(sampInfo0))
sampList = sampList[-1]

#- Format gene counts
for(i in sampList){
  tmp = read.table(file=paste("./RNAStarMulti/",i,".ReadsPerGene.out.tab",sep=""),sep="\t",header=F)
  tmp = tmp[-c(1:4),c(1,4)]
  colnames(tmp) = c('geneID',i)
  countData = merge(countData,tmp,by='geneID')
  rm(tmp)
}
sheepList = factor(sampInfo$sheep)
rownames(countData) = countData$geneID
countData$geneID = NULL
countData = countData[which(rowSums(countData)!=0),]

##- Get rid of 229_1 / outlier on PCA
countData = countData[,-1]
sampInfo = sampInfo0[-1,]

##--- dataframe creation
y0 = DGEList(counts = countData, 
             group=sampInfo$condition
)
dim(y0)
#[1] 18238    23
min(y0$samples$lib.size)
# 8466367

design <- model.matrix(~ condition, data = sampInfo)

##--- Filter genes with low counts
keep <- filterByExpr(y0, design = design, min.count = 5, min.total.count = 30)

y <- y0[keep,]
dim(y)
# [1] 12648     23

##-- dataframe creation = the one created for edgeR
dge <- y

##- Normalization
dge <- calcNormFactors(dge)
v <- voom(dge,design,plot=TRUE)

vfit <- lmFit(v,design)
vfit <- eBayes(vfit)

topVoom1 <- topTable(vfit,coef=ncol(design),number=Inf,sort.by="P")
sum(topVoom1$adj.P.Val<0.05) 
#[1] 112
summary(decideTests(vfit,adjust.method = "BH", p.value=0.05))
#        (Intercept) conditionS
# Down           985         89
# NotSig         583      12536
# Up           11080         23

vmlfc1 <- topTable(vfit,coef=ncol(design),number=Inf,sort.by="P",
                    adjust.method = "BH", p.value=0.05, lfc=1)
summary(decideTests(vfit,adjust.method = "BH", p.value=0.05, lfc=1))
#        (Intercept) conditionS
# Down           527         25
# NotSig        1566      12622
# Up           10555          1

vmlfc2 <- topTable(vfit,coef=ncol(design),number=Inf,sort.by="P",
                    adjust.method = "BH", p.value=0.05, lfc=2)
summary(decideTests(vfit,adjust.method = "BH", p.value=0.05, lfc=2))
#        (Intercept) conditionS
# Down            21          4
# NotSig        3031      12644
# Up            9596          0

selected.voom = factor(rownames(vmlfc1))
selected.voom
# [1] HCON_00191400 HCON_00182130 HCON_00193080 HCON_00007272 HCON_00132830 HCON_00121620 HCON_00192750
# [8] HCON_00193430 HCON_00029050 HCON_00050580 HCON_00090510 HCON_00191510 HCON_00115500 HCON_00007274
# [15] HCON_00053770 HCON_00106640 HCON_00134140 HCON_00191530 HCON_00090580 HCON_00191700 HCON_00087240
# [22] HCON_00035530 HCON_00157290 HCON_00048340 HCON_00017240 HCON_00074200

voom05 = topTable(vfit,coef=ncol(design),number=Inf,sort.by="P",
                  adjust.method = "BH", p.value=0.05)
dim(voom05)
#[1] 112  6
rm(y,y0)

#------=====#### SALMON counts analysis // VOOM #####=====------

##-- 2291_quant does not behave like other samples (38% variance): removed from further analyses
sampInfo = sampInfo0[-1,]
files <- file.path("./RNAsalm",paste0(sampInfo$sample,"_quant"), "quant.sf")
all(file.exists(files))
design <- model.matrix(~ condition, data = sampInfo)

##-- import Salmon files
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
names(txi)

y0 <- DGEList(txi$counts)
dim(y0)
# [1] 19384    23

##- filtering
keep <- filterByExpr(y0, design = design, min.count = 5, min.total.count = 15)
y <- y0[keep, ]
dim(y)
#[1] 10833    23

##- Normalization
y <- calcNormFactors(y)
v <- voom(y, design)
# v is now ready for lmFit() see limma User's Guide

###----- Analysis 
vfit.salm <- lmFit(v,design)
vfit.salm <- eBayes(vfit.salm)
topTable(vfit.salm,coef="conditionS",sort.by="P")
#                    logFC    AveExpr         t      P.Value  adj.P.Val          B
# HCON_00182130 -1.8208267  9.5644406 -5.824952 3.831747e-06 0.04150931  4.3404068
# HCON_00087240 -3.3700229  0.5623572 -5.387211 1.202655e-05 0.06514181  1.3889460
# HCON_00191700 -2.9444608  4.3946240 -5.027788 3.100383e-05 0.08093667  1.9251225
# HCON_00092800 -0.7697781 10.1673730 -5.008378 3.263476e-05 0.08093667  2.4042929
# HCON_00029050 -1.5257958  4.4822889 -4.925087 4.066904e-05 0.08093667  1.7611690
# HCON_00132830 -1.5547144  4.4894205 -4.888256 4.482784e-05 0.08093667  1.6859363
# HCON_00007272 -1.8122678  5.7787301 -4.662889 8.135002e-05 0.12589497  1.3926968
# HCON_00007274 -4.3978545  4.0314640 -4.409675 1.587954e-04 0.21502884  0.5600076
# HCON_00063030 -1.2800911  3.2073243 -4.329650 1.960878e-04 0.22010138  0.3129742
# HCON_00138210 -2.4904799  0.1479615 -4.316166 2.031768e-04 0.22010138 -0.3873761

topVoom.salm = topTable(vfit.salm,coef="conditionS",sort.by="P",number=Inf)
sum(topVoom.salm$adj.P.Val<0.05) 
# [1] 1
summary(vm.salm<-decideTests(vfit.salm,adjust.method = "BH", p.value=0.05))
#        (Intercept) conditionS
# Down             6          1
# NotSig        1024      10832
# Up            9803          0

vmlfc1.salm <- topTable(vfit.salm,coef=ncol(design),number=Inf,sort.by="P",
                    adjust.method = "BH", p.value=0.05, lfc=1)
summary(vm.salm<-decideTests(vfit.salm,adjust.method = "BH", p.value=0.05, lfc=1))
#        (Intercept) conditionS
# Down             4          1
# NotSig        1268      10832
# Up            9561          0

vmlfc2.salm <- topTable(vfit.salm,coef=ncol(design),number=Inf,sort.by="P",
                    adjust.method = "BH", p.value=0.05, lfc=2)
summary(vm2.salm <- decideTests(vfit.salm,adjust.method = "BH", p.value=0.05, lfc=2))
#        (Intercept) conditionS
# Down             0          0
# NotSig        2614      10833
# Up            8219          0

selected.voom.salm = factor(rownames(vmlfc1.salm))
selected.voom.salm
#[1] HCON_00182130

voom05.salm = topTable(vfit.salm,coef=ncol(design),number=Inf,sort.by="P",
                  adjust.method = "BH", p.value=0.05)
dim(voom05.salm)
#[1] 1  6

###-------====== Conclusions ======--------

###-- Summary
length(selected.salm)
#[1] 18
length(selected.star)
#[1] 38
length(selected.voom)
#[1] 26
length(selected.voom.salm)
#[1] 1

###-- Intersecting gene set Salmon DESeq2 / STAR DESeq2 / STAR VOOM / Salmon VOOM
df1 = data.frame(table(c(as.character(selected.voom),
                         as.character(selected.voom.salm),
                        as.character(selected.salm),
                        as.character(selected.star))))
dim(df1)
#[1] 48  2
df1$salm = match(df1$Var1,selected.salm)
df1$star = match(df1$Var1,selected.star)
df1$vostar = match(df1$Var1,selected.voom)
df1$vosalm = match(df1$Var1,selected.voom.salm)
df1$salm[is.na(df1$salm)] = 0
df1$star[is.na(df1$star)] = 0
df1$vostar[is.na(df1$vostar)] = 0
df1$vosalm[is.na(df1$vosalm)] = 0
df1$salm[df1$salm!=0] = 1
df1$star[df1$star!=0] = 1
df1$vostar[df1$vostar!=0] = 1
df1$vosalm[df1$vosalm!=0] = 1

##-- Conserved across at least 3 modalities
dfin = df1[df1$Freq >= 3,]
dfin$Var1
# [1] HCON_00007272 HCON_00007274 HCON_00029050 HCON_00050580 HCON_00087240 HCON_00090510 HCON_00090580
# [8] HCON_00132830 HCON_00182130 HCON_00191700

dfin
#             Var1 Freq salm star vostar vosalm
# 2  HCON_00007272    3    1    1      1      0
# 3  HCON_00007274    3    1    1      1      0
# 12 HCON_00029050    3    1    1      1      0
# 18 HCON_00050580    3    1    1      1      0
# 24 HCON_00087240    3    1    1      1      0
# 25 HCON_00090510    3    1    1      1      0
# 26 HCON_00090580    3    1    1      1      0
# 36 HCON_00132830    3    1    1      1      0
# 39 HCON_00182130    4    1    1      1      1
# 44 HCON_00191700    3    1    1      1      0

####--- Countdata for each gene
m2 = m1[m1$variable %in% dfin$Var1,]
head(m2)
# Group  samp      variable value
# 1     R 229_2 HCON_00007272  3807
# 2     R 229_3 HCON_00007272  3294
# 3     S 261_4 HCON_00007272   906
# 4     S 261_5 HCON_00007272  1139
# 5     S 261_6 HCON_00007272  1105
# 6     R 313_7 HCON_00007272 29490

aggregate(value ~ Group + variable, data= m2,FUN=mean)
#    Group      variable        value
# 1      R HCON_00007272  8557.818182
# 2      S HCON_00007272  1101.666667
# 3      R HCON_00007274 13519.818182
# 4      S HCON_00007274   235.916667
# 5      R HCON_00029050   933.090909
# 6      S HCON_00029050   376.500000
# 7      R HCON_00050580  6880.363636
# 8      S HCON_00050580  2181.166667
# 9      R HCON_00087240   107.909091
# 10     S HCON_00087240     4.916667
# 11     R HCON_00090510   109.818182
# 12     S HCON_00090510    24.916667
# 13     R HCON_00090580   114.636364
# 14     S HCON_00090580    40.166667
# 15     R HCON_00132830   464.363636
# 16     S HCON_00132830   148.750000
# 17     R HCON_00182130 50335.090909
# 18     S HCON_00182130 12979.083333
# 19     R HCON_00191700  2759.818182
# 20     S HCON_00191700   731.083333

ggplot(m2,aes(x = variable, y = value, col= Group)) +
  scale_y_log10() +
  geom_boxplot() + coord_flip()

####--- Summary stats for each gene

## VOOM-Star
vo = data.frame(vmlfc1[rownames(vmlfc1) %in% dfin$Var1,])
vo$gene = row.names(vo)
vo
#                   logFC    AveExpr         t      P.Value    adj.P.Val        B          gene
# HCON_00182130 -1.642327 10.3752420 -7.581354 7.488138e-08 0.0003495162 7.902324 HCON_00182130
# HCON_00007272 -2.295183  7.1470749 -7.207010 1.778404e-07 0.0004498651 7.348801 HCON_00007272
# HCON_00132830 -1.467092  3.7699369 -6.923035 3.469108e-07 0.0006488361 6.678441 HCON_00132830
# HCON_00029050 -1.173702  4.9729266 -6.353617 1.364601e-06 0.0013903982 5.422307 HCON_00029050
# HCON_00050580 -1.320190  7.6419420 -6.304340 1.539024e-06 0.0013903982 5.305556 HCON_00050580
# HCON_00090510 -1.988566  1.3995679 -5.949648 3.686258e-06 0.0029139872 3.955170 HCON_00090510
# HCON_00007274 -4.059951  5.5405333 -5.696109 6.934858e-06 0.0039869131 3.879871 HCON_00007274
# HCON_00090580 -1.391130  1.7739663 -4.851773 5.886439e-05 0.0147933038 1.780970 HCON_00090580
# HCON_00191700 -1.768043  5.9855089 -4.838977 6.082003e-05 0.0147933038 1.777432 HCON_00191700
# HCON_00087240 -3.773054 -0.1794151 -4.562935 1.231803e-04 0.0221518843 0.388617 HCON_00087240

summary(2^(abs(vo$logFC)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.256   2.658   3.264   5.589   4.673  16.679 

## VOOM-Salmon
vo.salm = data.frame(vmlfc1.salm[rownames(vmlfc1.salm) %in% dfin$Var1,])
vo.salm$gene = row.names(vo.salm)
vo.salm
#                   logFC  AveExpr         t      P.Value  adj.P.Val        B          gene
# HCON_00182130 -1.820827 9.564441 -5.824952 3.831747e-06 0.04150931 4.340407 HCON_00182130

## Salmon-DESeq2
sa = data.frame(sigFC1.salm[rownames(sigFC1.salm) %in% dfin$Var1,])
sa$gene = row.names(sa)
sa
#                 baseMean log2FoldChange     lfcSE      stat       pvalue         padj          gene
# HCON_00007274  566.16699      -5.246560 0.6844301 -7.665589 1.780121e-14 2.380556e-10 HCON_00007274
# HCON_00087240   19.87342      -3.664607 0.5940266 -6.169097 6.868128e-10 4.592374e-06 HCON_00087240
# HCON_00182130 4729.28008      -1.720016 0.3152093 -5.456742 4.849496e-08 2.161744e-04 HCON_00182130
# HCON_00132830  147.75536      -1.587756 0.2988219 -5.313386 1.076071e-07 3.597573e-04 HCON_00132830
# HCON_00050580 1725.99997      -1.711843 0.3579625 -4.782186 1.733994e-06 3.864783e-03 HCON_00050580
# HCON_00029050  144.15591      -1.418968 0.3073466 -4.616834 3.896381e-06 5.789590e-03 HCON_00029050
# HCON_00090510   10.94115      -1.989607 0.4742590 -4.195190 2.726428e-05 2.604323e-02 HCON_00090510
# HCON_00007272  368.40361      -1.532756 0.3704001 -4.138108 3.501814e-05 2.876348e-02 HCON_00007272
# HCON_00191700  235.25477      -2.356448 0.5669471 -4.156380 3.233295e-05 2.876348e-02 HCON_00191700
# HCON_00090580   19.92084      -1.984206 0.4936461 -4.019490 5.832425e-05 3.899851e-02 HCON_00090580

summary(sa$log2FoldChange)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -5.247  -2.265  -1.852  -2.321  -1.619  -1.419 

## STAR-DESeq2
st = data.frame(sigFC1[rownames(sigFC1) %in% m2$variable,])
st$gene = row.names(st)
st
#                  baseMean log2FoldChange     lfcSE      stat       pvalue         padj          gene
# HCON_00132830   270.79560      -1.569001 0.1905173 -8.235477 1.788432e-16 2.436917e-12 HCON_00132830
# HCON_00007272  4055.48471      -2.675834 0.3317620 -8.065524 7.292244e-16 3.312137e-12 HCON_00007272
# HCON_00182130 28002.74356      -1.766104 0.2212686 -7.981720 1.443083e-15 4.915863e-12 HCON_00182130
# HCON_00007274  5098.91744      -5.224254 0.6864959 -7.610029 2.740345e-14 6.223323e-11 HCON_00007274
# HCON_00087240    46.22323      -4.175536 0.5671727 -7.362019 1.811490e-13 3.085420e-10 HCON_00087240
# HCON_00090510    61.15105      -2.052257 0.2969646 -6.910780 4.819947e-12 6.567660e-09 HCON_00090510
# HCON_00050580  4003.95112      -1.492361 0.2203339 -6.773179 1.259826e-11 1.560580e-08 HCON_00050580
# HCON_00029050   587.10320      -1.173678 0.1877340 -6.251812 4.057176e-10 3.455193e-07 HCON_00029050
# HCON_00090580    69.82484      -1.497693 0.2639559 -5.674028 1.394784e-08 8.263187e-06 HCON_00090580
# HCON_00191700  1501.16021      -1.760377 0.3501690 -5.027220 4.976406e-07 1.832663e-04 HCON_00191700

summary(st$log2FoldChange)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -5.224  -2.520  -1.763  -2.339  -1.516  -1.174

## Correlation between DESEq2 STAR vs DESEq2 Salmon
sast = merge(sa,st,by='gene')
ggplot(sast,aes(x=log2FoldChange.x,y=log2FoldChange.y)) +
  geom_point() + 
  geom_abline(slope=1,intercept = 0,col='grey') +
  xlab('Salmon - log2FC') + 
  ylab('Star - log2FC') 

## Correlation between DESeq2/VOOM (STAR mapper)
stvo = merge(vo,st,by='gene')
ggplot(stvo,aes(x=logFC,y=log2FoldChange)) +
  geom_point() + 
  geom_abline(slope=1,intercept = 0,col='grey') +
  xlab('VOOM - log2FC') + 
  ylab('DESEq2 - log2FC') 

## Summary Fold Change based on STAR DESeq2
summary(2^(abs(st$log2FoldChange[match(st$gene,dfin$Var1)])))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.256   2.860   3.395   8.364   5.829  37.382  

###-- Summarize gene information from DESeq2-STAR analysis
m2$Sheep = factor(sapply(str_split(m2$samp,"_"),function(x)x[1]))
head(m2)
#   Group  samp      variable value Sheep
# 1     R 229_2 HCON_00007272  3807   229
# 2     R 229_3 HCON_00007272  3294   229
# 3     S 261_4 HCON_00007272   906   261
# 4     S 261_5 HCON_00007272  1139   261
# 5     S 261_6 HCON_00007272  1105   261
# 6     R 313_7 HCON_00007272 29490   313

library(scales)

pdf(file=paste0(draft,'./Supplementary_FigDEgenes.pdf'))
p1 = ggplot(m2,aes(x=variable,y=value, shape=Group, color=Sheep, group=Group)) +
  geom_point(size = 4,alpha = .3,position=position_dodge(width=.8)) +
  annotation_logticks(scaled = TRUE,sides="l") +
  scale_y_log10(breaks = c(0.000001, 10^(-2:6)), 
                labels = c(0, math_format()(-2:6))) +
  ylab('Gene-wise transcript count') + xlab('Gene name') +
  theme_bw() + 
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, vjust=0.5),legend.position = 'bottom')
dev.off()

####--- Correlation with FEC and Ht
#setwd('~/Documents/INRA/GEMANEMA/Analysis/')
a = read.csv(file='../parasito_infestation1et2.csv',header=T,sep=",")
#f = read.csv(file="../indexASC&SNP_180G1.csv",sep=";",header=T)

a$Sheep = sapply(str_split(a$IPG,'20000152'), function(x) x[2])

#f$Sheep = sapply(str_split(f$IPG,'20000152'), function(x) x[2])
a$Sheep[a$Sheep==232] = 432 ## rename sheep
#f$Sheep[f$Sheep==232] = 432 ## rename sheep
pheno = merge(m2,a[,c('opg302','ht302','Sheep')],by='Sheep')
#pheno = merge(pheno, f[,c('sol_opg' ,'indexSNP','Sheep')],by='Sheep')

##-- Relationship between transcript count and FEC
colnames(pheno)[4] ='Gene'
p_rnafec = ggplot(pheno,aes(x = opg302,y = value,col=Gene,fill=Gene)) + 
  facet_wrap(~ Gene,ncol=3) +
  annotation_logticks(scaled = TRUE,sides="l") +
  scale_y_log10(breaks = c(0.000001, 10^(-2:5)), 
                labels = c(0, math_format()(-2:5))) +
  scale_x_log10(breaks = c(0.000001, 10^(-2:4)), 
                labels = c(0, math_format()(-2:4))) +
  geom_point(shape=pheno$Group,size=3) + geom_smooth(method = 'lm',alpha=.3,lwd=.5) +
  theme(legend.position = 'none',text = element_text(size = 14)) + 
  xlab('Faecal Egg Count (eggs/g)') + ylab('Gene-wise transcript count')

n = 0
mat = matrix(0,dim(dfin)[1],3)
for(g in unique(pheno$Gene)){
  n = n+1
  mat[n,1] = g
  pears = rcorr(log(pheno$opg302[pheno$Gene==g]+1),log(pheno$value[pheno$Gene==g]+1))
  mat[n,2] = pears$r[1,2]
  mat[n,3] = pears$P[1,2]
}

mat = data.frame(mat)
colnames(mat) = c('gene','Correlation','P')
mat
#             gene        Correlation                    P
# 1  HCON_00007272 -0.876552166809371 4.15159244759877e-08
# 2  HCON_00132830 -0.787005942792415 8.39257679197125e-06
# 3  HCON_00050580 -0.797399399918503 5.21460649416383e-06
# 4  HCON_00090580 -0.738821231744029 5.66319265344895e-05
# 5  HCON_00029050 -0.751864666670976 3.52336519142682e-05
# 6  HCON_00191700  -0.74449960256976 4.62193166901326e-05
# 7  HCON_00090510 -0.847828450221167 3.27101198349666e-07
# 8  HCON_00007274  -0.86589320180777 9.43269511388678e-08
# 9  HCON_00087240 -0.876317844283727 4.23054666853773e-08
# 10 HCON_00182130 -0.859005491216781 1.54614228531358e-07

summary(as.numeric(as.character(mat$Correlation)))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.8766 -0.8642 -0.8226 -0.8145 -0.7607 -0.7388 

###--- Figure output
pdf(file=paste0(draft,'supplementary_Figure5.pdf'),width=10,height=12)
print(prna_fec) 
dev.off()

# ##-- Relationship between transcript count and Hct
# p3 = ggplot(pheno,aes(x = ht302,y = value,col=Gene,fill=Gene)) + 
#   facet_wrap(~ Gene ) +
#   annotation_logticks(scaled = TRUE,sides="l") +
#   scale_y_log10(breaks = c(0.000001, 10^(-2:5)), 
#                 labels = c(0, math_format()(-2:5))) +
#   geom_point(shape=pheno$Group,size=3) + geom_smooth(method = 'lm',alpha=.3,lwd=.5) +
#   theme(legend.position = 'none') + xlab('Haematocrit (%)') + ylab('Gene-wise transcript count')
# 
# n = 0
# mat2 = matrix(0,dim(dfin)[1],3)
# for(g in unique(pheno$Gene)){
#   n = n+1
#   mat2[n,1] = g
#   pears = rcorr(log(pheno$ht302[pheno$Gene==g]+1),log(pheno$value[pheno$Gene==g]+1))
#   mat2[n,2] = pears$r[1,2]
#   mat2[n,3] = pears$P[1,2]
# }
# 
# mat2 = data.frame(mat2)
# colnames(mat2) = c('gene','Correlation','P')
# mat2
# #             gene       Correlation                    P
# # 1  HCON_00007272 0.579127516941794  0.00378381803362826
# # 2  HCON_00132830 0.617893097098845  0.00167850232588362
# # 3  HCON_00050580 0.560234012098864  0.00543089716332257
# # 4  HCON_00090580 0.506514923922898   0.0136491633322109
# # 5  HCON_00029050 0.466320833796822    0.024902193248749
# # 6  HCON_00191700 0.591276558055554  0.00296482351977079
# # 7  HCON_00090510  0.66707064133416 0.000507619144310301
# # 8  HCON_00007274 0.587182708054064  0.00322218829661969
# # 9  HCON_00087240 0.657891670053722 0.000644811564712633
# # 10 HCON_00182130 0.598101501907202  0.00257434451312455
# 
# summary(as.numeric(as.character(mat2$Correlation)))
# #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.4663  0.5650  0.5892  0.5832  0.6129  0.6671 

# pdf(file=paste0(draft,'supplementary_Figure6.pdf'),width=10,height=12)
# print(p3) 
# dev.off()

####--- Retrieve gene expression associated with regulators of collagen exp:
# dpy-7:HCON_00190185 
# skn-1:HCON_00096900
# elt-3:HCON_00035090 / HCON_00084430 / HCON_00188830

cutreg = c('HCON_00190185','HCON_00096900','HCON_00035090','HCON_00084430','HCON_00188830')
match(cutreg,rownames(voom05))


####--- qRT-PCR 
setwd('../RTPCR/Validation/')

##--- Melting curves
glist1 = c('ama','far','gpd','HCON_00007272','HCON_00007272_1','HCON_00007274','HCON_00007274_1',
          'HCON_00029050','HCON_00029050_1','HCON_00050580','HCON_00050580_1','empty')
glist2 = c('HCON_00087240','HCON_00087240_1','HCON_00090510','HCON_00090510_1',
           'HCON_00090580','HCON_00090580_1','HCON_00132830','HCON_00132830_1',
           'HCON_00182130_1','HCON_00182130','HCON_00191700','HCON_00191700_1')
grad = data.frame(let = c('C', 'D', 'E'), temp = c(65.3,63.4,61))

## Plate 1
mdc1 = read.csv(file ='./MPN_2019-12-16_GEMANEMA_Gradient1 -  Melt Curve Derivative Results_SYBR.csv',header = T,sep=';')
ac1 = read.csv(file ='./MPN_2019-12-16_GEMANEMA_Gradient1 -  Melt Curve RFU Results_SYBR.csv',header = T,sep=';')
acy1 = read.csv(file = './MPN_2019-12-16_GEMANEMA_Gradient1 -  Quantification Amplification Results_SYBR.csv',header = T,sep=';')

## Plate 2
mdc2 = read.csv(file ='./MPN_2019-12-16_Gradient2 -  Melt Curve Derivative Results_SYBR.csv',header = T,sep=';')
ac2 = read.csv(file ='./MPN_2019-12-16_Gradient2 -  Melt Curve RFU Results_SYBR.csv',header = T,sep=';')
acy2 = read.csv(file = './MPN_2019-12-16_GEMANEMA_Gradient1 -  Quantification Amplification Results_SYBR.csv',header = T,sep=';')

convert_mdc<- function(mdc,gene.list){
  mdc = reshape2::melt(mdc,1)
  ## Keep wells of interest C D & E
  mdc = mdc[c(grep('C',mdc$variable),grep('D',mdc$variable),grep('E',mdc$variable)),]
  mdc$variable = factor(mdc$variable)
  mdc$gene = rep(gene.list,each = 83, times = 3)
  return(mdc)
}
convert_acy<- function(acy,gene.list){
  acy = reshape2::melt(acy,1)
  ## Keep wells of interest C D & E
  acy = acy[c(grep('C',acy$variable),grep('D',acy$variable),grep('E',acy$variable)),]
  acy$variable = factor(acy$variable)
  acy$gene = rep(gene.list,each = 40, times = 3)
  acy$temp = grad$temp[match(substr(acy$variable,1,1),grad$let)]
  return(acy)
}

mdc1.convert = convert_mdc(mdc1,glist1)
ac1.convert =  convert_mdc(ac1,glist1)
acy1.convert = convert_acy(acy1,glist1)

mdc2.convert = convert_mdc(mdc2,glist2)
ac2.convert =  convert_mdc(ac2,glist2)
acy2.convert = convert_acy(acy2,glist2)

## Plot 1st regression of relative fluorescence vs. temperature
p1 = ggplot(mdc1.convert,aes(x = Temperature, y = value,gp = variable,col=gene)) +
  geom_line() + facet_wrap(~ gene) +
  ylab('-dRFU/dT') + xlab('Temperature')  + 
  theme(legend.position = 'none',text = element_text(size = 8))

## Plot normalized fluorescence vs. temperature
p2 = ggplot(ac1.convert,aes(x = Temperature, y = value,gp = variable,col=gene)) +
  geom_line() + facet_wrap(~ gene) +
  ylab('Normalized RFU') + xlab('Temperature') + 
  theme(legend.position = 'none',text = element_text(size = 8))

## Plot 1st regression of relative fluorescence vs. temperature
p3 = ggplot(mdc2.convert,aes(x = Temperature, y = value,gp = variable,col=gene)) +
  geom_line() + facet_wrap(~ gene) +
  ylab('-dRFU/dT') + xlab('Temperature')  + 
  theme(legend.position = 'none',text = element_text(size = 8))

## Plot normalized fluorescence vs. temperature
p4 = ggplot(ac2.convert,aes(x = Temperature, y = value,gp = variable,col=gene)) +
  geom_line() + facet_wrap(~ gene) +
  ylab('Normalized RFU') + xlab('Temperature') + 
  theme(legend.position = 'none',text = element_text(size = 8))

multiplot(p1,p2,p3,p4,cols=2)

acy.convert = rbind(acy1.convert,acy2.convert)

pdf(file = './BestTemp.pdf')
ggplot(acy.convert,aes(x = Cycle, y = value,gp = variable,col=factor(temp))) +
  geom_line() + facet_wrap(~ gene) +
  ylab('Normalized RFU') + xlab('Temperature') + 
  theme(legend.position = 'bottom',text = element_text(size = 8))
dev.off()

##--- Efficacy curves

## create metadata
cc = c(1/10,1/50,1/250,1/1250)
glist1.1 = c('ama','far','gpd','HCON_00007272')
glist1.2 = c('HCON_00007272_1','HCON_00007274','HCON_00029050_1','HCON_00050580')

glist2.1 = c('HCON_00050580_1','HCON_00087240_1','HCON_00090510','empty')
glist2.2 = c('HCON_00090510_1','HCON_00090580','HCON_00090580_1','HCON_00132830')

glist3.1 = c('HCON_00132830_1','HCON_00182130_1','HCON_00182130','empty')
glist3.2 = c('HCON_00191700','HCON_00191700_1','empty','empty')

let = c('A','B','C','D','E','F','G','H')
v = NULL
n=0
for(i in let){
  for(j in seq(1:12)){
    n=n+1
    v[n] = paste0(i,j)
  }
}

metaqrt = data.frame(plate = rep(c(1,2,3),each = 96),
                     conc = rep(rep(cc, each = 12, times = 2),3),
                     gene = c(rep(glist1.1,each = 3, times = 4),
                              rep(glist1.2,each = 3, times = 4),
                              rep(glist2.1,each = 3, times = 4),
                              rep(glist2.2,each = 3, times = 4),
                              rep(glist3.1,each = 3, times = 4),
                              rep(glist3.2,each = 3, times = 4)),
                     well = rep(v,3))

## Plate 1
ac1 = read.csv(file ='./MPN_2019-12-17_GEMANEMA_Efficacy1 -  Quantification Amplification Results_SYBR.csv',header = T,sep=';')
cq1.0 = read.csv(file ='./MPN_2019-12-17_GEMANEMA_Efficacy1 -  Quantification Cq Results_0.csv',header = T,sep=';')

## Plate 2
ac2 = read.csv(file ='./MPN_2019-12-17_GEMANEMA_Efficacy2 -  Quantification Amplification Results_SYBR.csv',header = T,sep=';')
cq2.0 = read.csv(file ='./MPN_2019-12-17_GEMANEMA_Efficacy2 -  Quantification Cq Results_0.csv',header = T,sep=';')

## Plate 3
ac3 = read.csv(file ='./MPN_2019-12-17_GEMANEMA_Efficacy3 -  Quantification Amplification Results_SYBR.csv',header = T,sep=';')
cq3.0 = read.csv(file ='./MPN_2019-12-17_GEMANEMA_Efficacy3 -  Quantification Cq Results_0.csv',header = T,sep=';')

## Retrieve amplification curves
convert_ac<- function(ac,meta,plt){
  ac = reshape2::melt(ac1,1)
  colnames(ac)[2] = 'well'
  ## Add metadata
  ac = merge(ac,meta[meta$plate == plt,],by='well')
  return(ac)
  rm(ac)
}

ac1.convert = convert_ac(ac1,metaqrt,1)
ac2.convert = convert_ac(ac2,metaqrt,2)
ac3.convert = convert_ac(ac3,metaqrt,3)

dim(ac1.convert)
#[1] 3840    6
dim(ac2.convert)
#[1] 3840    6
dim(ac3.convert)
#[1] 3840    6

## Retrieve Ct data
cq1 = data.frame(Cq = cq1.0$Cq, well = v)
cq1 = merge(cq1, metaqrt[metaqrt$plate==1,], by = 'well')

cq2 = data.frame(Cq = cq2.0$Cq, well = v)
cq2 = merge(cq2, metaqrt[metaqrt$plate==2,], by = 'well')

cq3 = data.frame(Cq = cq3.0$Cq, well = v)
cq3 = merge(cq3, metaqrt[metaqrt$plate==3,], by = 'well')

####--- Plot data
ac.convert = rbind(ac1.convert, ac2.convert, ac3.convert)
cq = rbind(cq1, cq2, cq3)
cq$Cq = as.numeric(as.character(cq$Cq))

pdf(file = './EfficacyCurve.pdf')
ggplot(cq[cq$gene!='empty',], aes(x = log10(conc), y = Cq)) +
  geom_smooth(method = 'lm',se = T, col = 'blue',fill = 'grey', lwd = 0.8) +
  geom_point(alpha = .3, col = 'red') + 
  facet_wrap(~ gene) + #, scales = 'free') + 
  #scale_y_log10() + 
  #scale_x_log10() + 
  xlab('cDNA dilution') +
  theme(text = element_text(size = 8))
dev.off()

pdf(file = './BestConc.pdf')
ggplot(ac.convert[ac.convert$gene!='empty',], aes(x = Cycle, y = value, gp = well, col = factor(conc))) +
  geom_line() + 
  facet_wrap(~ gene, scales = 'free') + 
  scale_x_log10() + xlab('Cycle') + ylab('RFU') +
  #theme(legend.position = 'bottom')
  theme(text = element_text(size = 8), legend.position = 'bottom')
dev.off()

####--- Compute efficacies
#Slopes between -3.1 and -3.6 giving reaction efficiencies between 90 and 110% are typically acceptable. 
eff = NULL
cq = cq[cq$gene!='empty',]
for(g in unique(cq$gene)){
  s = summary(lm(Cq ~ log10(conc), data = cq[cq$gene ==g,]))
  tmp = data.frame(gene = g,
                   slope = s$coefficients[2,1],
                   slope.se = s$coefficients[2,2],
                   E = 10^(-1/s$coefficients[2,1])-1,
                   interc = s$coefficients[1,1],
                   R2 = s$adj.r.squared)
  eff = rbind(eff,tmp)
}

eff[order(as.character(eff$gene)),]
#               gene slope slope.se     E interc    R2
# 1              ama -3.12   0.1103 1.092   23.9 0.986
# 3              far -2.91   0.0853 1.207   16.9 0.991
# 4              gpd -3.01   0.0714 1.149   16.1 0.994
# 2    HCON_00007272 -3.19   0.1394 1.060   20.4 0.979
# 5  HCON_00007272_1 -3.23   0.1477 1.042   21.1 0.977
# 7    HCON_00007274 -3.30   0.0823 1.011   20.7 0.993
# 8  HCON_00029050_1 -3.23   0.0777 1.040   19.3 0.994
# 6    HCON_00050580 -3.16   0.1082 1.072   20.3 0.987
# 9  HCON_00050580_1 -2.93   0.1710 1.195   19.9 0.964
# 10 HCON_00087240_1 -3.44   0.4310 0.953   24.7 0.851
# 11   HCON_00090510 -1.57   0.1542 3.322   27.2 0.904
# 12 HCON_00090510_1 -1.57   0.1294 3.327   25.7 0.930
# 14   HCON_00090580 -3.02   0.1007 1.144   22.2 0.988
# 15 HCON_00090580_1 -3.28   0.1246 1.016   21.2 0.984
# 13   HCON_00132830 -3.17   0.1270 1.066   21.5 0.983
# 16 HCON_00132830_1 -3.41   0.1337 0.966   20.6 0.985
# 18   HCON_00182130 -3.14   0.0774 1.084   18.3 0.993
# 17 HCON_00182130_1 -3.20   0.0688 1.055   18.6 0.995
# 19   HCON_00191700 -3.33   0.0941 0.998   20.5 0.991
# 20 HCON_00191700_1 -3.43   0.1698 0.955   20.6 0.974

###--- Quantification
diff1 = read.csv(file = '../Differential/MPN_2019-12-17_GEMANEMA_Differential1 -  Quantification Cq Results_0.csv', header = T,sep=';')
diff2 = read.csv(file = '../Differential/MPN_2019-12-17_GEMANEMA_Differential2 -  Quantification Cq Results_0.csv', header = T,sep=';')
diff3 = read.csv(file = '../Differential/MPN_2019-12-17_GEMANEMA_Differential3 -  Quantification Cq Results_0.csv', header = T,sep=';')
## Note: plate 4 involves a shift in columns (column 4 is empty but for ama; control wells were A12-F12)
diff4 = read.csv(file = '../Differential/MPN_2020-02-04_GEMANEMA_Differential1_PlateCorrected -  Quantification Cq Results_0.csv', header = T,sep=';')
## Note: plate 5: ama was loaded on A1-3 and B1-3 / far was on E1-3 and F1-3 in lieu; values replaced in the correct spots manually
## Note2: worms from the S sheep were loaded first
diff5 = read.csv(file = '../Differential/MPN_2020-02-04_GEMANEMA_Differential2_PlateCorrected -  Quantification Cq Results_0.csv', header = T,sep=';')
## Note: plate 6: Blank for R sample is duplicated A and B10-12
diff6 = read.csv(file = '../Differential/MPN_2020-02-04_GEMANEMA_Differential3 -  Quantification Cq Results_0.csv', header = T,sep=';')

diff1$Well = v ## convert A01 to A1
diff2$Well = v
diff3$Well = v
diff4$Well = v 
diff5$Well = v
diff6$Well = v

## Remove control and empty wells
emptyWell = c(grep(10,diff1$Well),
              grep(11,diff1$Well),
              grep(12,diff1$Well))

diff1 = diff1[which(!(seq(1:dim(diff1)[1]) %in% emptyWell)),c(1,7)]
diff2 = diff2[which(!(seq(1:dim(diff2)[1]) %in% emptyWell)),c(1,7)]
diff3 = diff3[which(!(seq(1:dim(diff3)[1]) %in% emptyWell)),c(1,7)]
diff4 = diff4[which(!(seq(1:dim(diff4)[1]) %in% emptyWell)),c(1,7)]
diff5 = diff5[which(!(seq(1:dim(diff5)[1]) %in% emptyWell)),c(1,7)]
diff6 = diff6[which(!(seq(1:dim(diff6)[1]) %in% emptyWell)),c(1,7)]

## Assign plate number
diff1$Plate = 1
diff2$Plate = 2
diff3$Plate = 3
diff4$Plate = 4
diff5$Plate = 5
diff6$Plate = 6

## Gene list for plates 1 to 3
glist = c('ama','far','gpd','HCON_00007272','HCON_00007274','HCON_00029050_1','HCON_00050580',
           'HCON_00087240_1','HCON_00090580_1','HCON_00132830','HCON_00182130','HCON_00191700')
hkg.list = c('ama','far','gpd')
samp = c('R1','S1','R2','S2','R3','S3','R4','S4','S5','R5','R6','S6')
let = c('A','B','C','D','E','F','G','H')
v2 = NULL
n=0
for(i in seq(1:9)){
  for(j in let){
    n=n+1
    v2[n] = paste0(j,i)
  }
}
## Metadata 
metadif = data.frame(Plate = rep(c(1,2,3,4,5,6),each = 72),
                     Well = rep(v2,6),
                     Gene = rep(c(rep(glist[1:8],3), rep(c(glist[9:12],glist[1:4]),3),rep(glist[5:12],3)),6),
                     Sample = c(rep(samp[1],24),rep(c(samp[1],samp[2]),each =4,times=3),rep(samp[2],24),
                                rep(samp[3],24),rep(c(samp[3],samp[4]),each =4,times =3),rep(samp[4],24),
                                rep(samp[5],24),rep(c(samp[5],samp[6]),each =4,times =3),rep(samp[6],24),
                                rep(samp[7],24),rep(c(samp[7],samp[8]),each =4,times=3),rep(samp[8],24),
                                rep(samp[9],24),rep(c(samp[9],samp[10]),each =4,times =3),rep(samp[10],24),
                                rep(samp[11],24),rep(c(samp[11],samp[12]),each =4,times =3),rep(samp[12],24))
                                )

## Merge
difftot = rbind(diff1,diff2,diff3,diff4,diff5,diff6)

diff = merge (difftot,metadif, by = c('Plate','Well'))

## S1 appears to be an outlier
aggregate(Cq ~ Sample, FUN = mean,data=diff)
#    Sample   Cq
# 1      R1 29.7
# 2      R2 30.1
# 3      R3 29.5
# 4      R4 29.6
# 5      R5 29.3
# 6      R6 28.8
# 7      S1 32.2
# 8      S2 30.6
# 9      S3 30.1
# 10     S4 29.6
# 11     S5 30.2
# 12     S6 28.9

## Compute reference gene stability according to Vandersompele 2002 (M factor)
M = NULL
for(j in hkg.list){
  V = NULL
  n=0
  for(k in hkg.list){
    if(j != k){ #} & !is.na(diff$Cq[diff$Gene==j]) & !is.na(diff$Cq[diff$Gene==k])){
      n = n+1
      nona = which(!is.na(diff$Cq[diff$Gene==j]) &!is.na(diff$Cq[diff$Gene==k]))
      A = log(diff$Cq[diff$Gene==j][nona]/diff$Cq[diff$Gene==k][nona])
      V[n] = sd(A)
    }
  }
  M[j] = sum(V)
  rm(A,V)
}
M
#    ama    far    gpd 
# 0.0537 0.0407 0.0431 

## Compute dCt: ama has higher expression level
hkg = aggregate(Cq ~ Gene + Sample, FUN = mean,data = diff[diff$Gene %in% hkg.list,])
hkg
#    Gene Sample   Cq
# 1   ama     R1 33.5
# 2   far     R1 26.5
# 3   gpd     R1 26.2
# 4   ama     R2 33.7
# 5   far     R2 26.8
# 6   gpd     R2 26.2
# 7   ama     R3 33.7
# 8   far     R3 26.3
# 9   gpd     R3 25.5
# 10  ama     S1 32.6
# 11  far     S1 29.3
# 12  gpd     S1 29.1
# 13  ama     S2 33.9
# 14  far     S2 26.9
# 15  gpd     S2 26.2
# 16  ama     S3 33.2
# 17  far     S3 26.1
# 18  gpd     S3 25.7

## Compute control avg Ct based on 3 hkg
hkg.m = aggregate(Cq ~ Sample, FUN = mean,data = hkg) ##ama does not behave as other genes
diff$hkg = hkg.m$Cq[match(diff$Sample,hkg.m$Sample)]
diff$Group = factor(substr(diff$Sample,1,1))

## Gene expression modeling
diff$Group = factor(diff$Group,levels = c('S','R'))
mCq = nlme::lme(Cq ~ Gene*Group + hkg,
          random =~ 1|Sample,
          data = na.omit(diff))
#plot(mCq)
s = summary(mCq)
s
# Linear mixed-effects model fit by REML
# Data: na.omit(diff) 
# AIC  BIC logLik
# 1204 1311   -575
# 
# Random effects:
#   Formula: ~1 | Sample
# (Intercept) Residual
# StdDev:       0.162    0.942
# 
# Fixed effects: Cq ~ Gene * Group + hkg 
#                            Value Std.Error  DF t-value p-value
# (Intercept)                11.11     1.813 384    6.13  0.0000
# Genefar                    -6.34     0.331 384  -19.16  0.0000
# Genegpd                    -6.93     0.335 384  -20.65  0.0000
# GeneHCON_00007272          -1.37     0.331 384   -4.15  0.0000
# GeneHCON_00007274          -3.64     0.340 384  -10.69  0.0000
# GeneHCON_00029050_1        -5.71     0.331 384  -17.27  0.0000
# GeneHCON_00050580          -2.46     0.335 384   -7.35  0.0000
# GeneHCON_00087240_1         0.73     0.345 384    2.11  0.0359
# GeneHCON_00090580_1        -4.87     0.331 384  -14.74  0.0000
# GeneHCON_00132830          -5.50     0.331 384  -16.62  0.0000
# GeneHCON_00182130          -3.43     0.335 384  -10.26  0.0000
# GeneHCON_00191700          -1.39     0.335 384   -4.16  0.0000
# GroupR                      0.22     0.348   9    0.62  0.5489
# hkg                         0.78     0.062   9   12.49  0.0000
# Genefar:GroupR             -0.85     0.459 384   -1.85  0.0651
# Genegpd:GroupR             -0.83     0.466 384   -1.77  0.0773
# GeneHCON_00007272:GroupR   -2.13     0.459 384   -4.63  0.0000
# GeneHCON_00007274:GroupR   -0.36     0.466 384   -0.78  0.4348
# GeneHCON_00029050_1:GroupR -0.09     0.459 384   -0.20  0.8401
# GeneHCON_00050580:GroupR   -1.52     0.462 384   -3.28  0.0011
# GeneHCON_00087240_1:GroupR  0.54     0.470 384    1.16  0.2478
# GeneHCON_00090580_1:GroupR  0.01     0.459 384    0.02  0.9868
# GeneHCON_00132830:GroupR    0.14     0.459 384    0.30  0.7668
# GeneHCON_00182130:GroupR   -2.41     0.462 384   -5.20  0.0000
# GeneHCON_00191700:GroupR   -1.63     0.462 384   -3.53  0.0005

r = s$tTable[which(s$tTable[,5]<0.05),]
r[grep(':',names(s$tTable[which(s$tTable[,5]<0.05),5])),]
#                          Value Std.Error  DF t-value  p-value
# GeneHCON_00007272:GroupR -2.13     0.459 384   -4.63 5.01e-06
# GeneHCON_00050580:GroupR -1.52     0.462 384   -3.28 1.13e-03
# GeneHCON_00182130:GroupR -2.41     0.462 384   -5.20 3.17e-07
# GeneHCON_00191700:GroupR -1.63     0.462 384   -3.53 4.72e-04

ggplot(diff,aes(x = Gene, y = Cq, col = Sample)) + 
  geom_point()

##---- Compute ddCt between R & S groups by gene
## relative to S group
## ∆∆Ct = ∆Ct (Sample) – ∆Ct (Control average)
## 1. Compute control dCt based on 3 hkg
require(dplyr)

diff.dCt = diff %>% 
  group_by(Sample,Gene) %>%
  mutate(dCt = Cq - hkg)

dCtS = aggregate(dCt ~ Gene, FUN = mean,data = diff.dCt[diff.dCt$Group=='S',]) ## Compute S group expression level for every gene
df = aggregate(dCt ~ Sample + Gene,FUN = mean,data = diff.dCt)  ## 
df$dCtS = dCtS$dCt[match(df$Gene,dCtS$Gene)]

df$ddCt = df$dCt - df$dCtS

## 2.Fold gene estimation = 2^-(∆∆Ct)
df$fc = 2^-(df$ddCt)
df$Group = factor(substr(df$Sample,1,1))

fc.tab = aggregate(fc ~ Group + Gene, data = df, FUN = mean)
fc.tab$std = aggregate(fc ~ Group + Gene, data = df, FUN = sd)[,3]
fc.tab[fc.tab$Group == 'R',]
#    Group            Gene    fc   std
# 1      R             ama 0.920 0.258
# 3      R             far 1.483 0.346
# 5      R             gpd 1.404 0.297
# 7      R   HCON_00007272 3.586 0.729
# 9      R   HCON_00007274 1.260 1.097
# 11     R HCON_00029050_1 0.867 0.131
# 13     R   HCON_00050580 2.455 0.702
# 15     R HCON_00087240_1 0.629 0.303
# 17     R HCON_00090580_1 0.875 0.374
# 19     R   HCON_00132830 0.750 0.178
# 21     R   HCON_00182130 4.513 1.267
# 23     R   HCON_00191700 2.711 1.070

## Plot
require(ggrepel)

ggplot(df,aes(x = Gene,y = fc,col = Group,label = Sample)) + 
  geom_point(size = 3,alpha = .4,position = position_dodge(0.5)) + 
  geom_text() + 
  scale_color_manual(values = colsRS) +
  scale_y_log10() + 
  geom_hline(yintercept = 1) +
  ylab('Fold gene relative to S worms') + xlab('Gene') +
  theme(axis.text.x = element_text(vjust = 0.6, angle = 90),
        text = element_text(size = 14))

res=NULL
for(i in 1:length(levels(df$Gene))){
  t = wilcox.test(fc ~ Group,data=df[df$Gene==levels(df$Gene)[i],])$p.value
  tmp = data.frame(p=t,g=levels(df$Gene)[i])
  res = rbind(res,tmp)
}
res[res$p<0.05/9,]
#          p             g
# 4  0.00216 HCON_00007272
# 7  0.00216 HCON_00050580
# 11 0.00216 HCON_00182130
# 12 0.00216 HCON_00191700

####---- S1 has an outlier behaviour: removed from data -- this does not alter the results but normalization is better
diff2 = diff[diff$Sample!='S1',]
hkg2 = aggregate(Cq ~ Gene + Sample, FUN = mean,data = diff2[diff2$Gene %in% hkg.list,])

hkg.m2 = aggregate(Cq ~ Sample, FUN = mean,data = hkg2)
diff2$hkg = hkg.m2$Cq[match(diff2$Sample,hkg.m2$Sample)]
diff2$Group = factor(substr(diff2$Sample,1,1))
diff2$Group = factor(diff2$Group,levels = c('S','R'))

## Linear model (underpowered )
mCq2 = nlme::lme(Cq ~ Gene*Group + hkg,
                 random =~ 1|Sample,
                 data = na.omit(diff2))
s2 = summary(mCq2)
r2 = s2$tTable[which(s2$tTable[,5]<0.05),]
r2[grep(':',names(s2$tTable[which(s2$tTable[,5]<0.05),5])),]
#                          Value Std.Error  DF t-value  p-value
# GeneHCON_00007272:GroupR -2.09     0.355 356   -5.90 8.48e-09
# GeneHCON_00050580:GroupR -1.41     0.355 356   -3.97 8.69e-05
# GeneHCON_00182130:GroupR -1.94     0.355 356   -5.47 8.46e-08
# GeneHCON_00191700:GroupR -1.54     0.355 356   -4.33 1.94e-05

## Gene expression modeling
##---- Compute ddCt between R & S groups by gene
## relative to S group
## ∆∆Ct = ∆Ct (Sample) – ∆Ct (Control average)
## 1. Compute control dCt based on 3 hkg
require(dplyr)

diff2.dCt = diff2 %>%
  group_by(Sample,Gene) %>%
  mutate(dCt = Cq - hkg)

dCtS2 = aggregate(dCt ~ Gene, FUN = mean,data = diff2.dCt[diff2.dCt$Group=='S',]) ## Compute S group expression level for every gene
df2 = aggregate(dCt ~ Sample + Gene,FUN = mean,data = diff2.dCt)  ##
df2$dCtS = dCtS2$dCt[match(df2$Gene,dCtS2$Gene)]

df2$ddCt2 = df2$dCt - df2$dCtS

## 2.Fold gene estimation = 2^-(∆∆Ct)
df2$fc = 2^-(df2$ddCt)
df2$Group = factor(substr(df2$Sample,1,1))

fc.tab = aggregate(fc ~ Group + Gene, data = df2, FUN = mean)
fc.tab$std = aggregate(fc ~ Group + Gene, data = df2, FUN = sd)[,3]
fc.tab[fc.tab$Group == 'R',]
#    Group            Gene    fc   std
# 1      R             ama 0.920 0.258
# 3      R             far 1.099 0.257
# 5      R             gpd 1.072 0.227
# 7      R   HCON_00007272 3.848 0.783 *
# 9      R   HCON_00007274 1.771 1.541
# 11     R HCON_00029050_1 1.292 0.196
# 13     R   HCON_00050580 2.436 0.696 *
# 15     R HCON_00087240_1 0.893 0.430
# 17     R HCON_00090580_1 1.119 0.478
# 19     R   HCON_00132830 1.120 0.265
# 21     R   HCON_00182130 3.497 0.982 *
# 23     R   HCON_00191700 2.715 1.072 *

df2$Gene=substr(df2$Gene,1,13)
ggplot(df2[!(df2$Gene %in% c('ama','far','gpd')),],
       aes(x = Gene,y = fc,col = Group,label = Sample)) +
  geom_point(size = 4,alpha = .4,position = position_dodge(.2)) +
  #geom_violin() +
  #geom_bar(stat) +
  #geom_text_repel() +
  scale_color_manual(values = colsRS) +
  #scale_fill_manual(values = colsRS) +
  scale_y_log10() + geom_hline(yintercept = 1) +
  theme_classic() +
  ylab('Fold gene relative to S worms') + xlab('Gene') +
  theme(axis.text.x = element_text(vjust = 0.6, angle = 90),
        text = element_text(size = 14))
dev.off()

##### Plot Figure 3

sig = c('HCON_00007272','HCON_00050580','HCON_00182130','HCON_00191700')

p_rnafec = ggplot(pheno[pheno$Gene %in% sig,],
                  aes(x = opg302,y = value,col=Group,fill=Group,group = Gene)) + 
  facet_wrap(~ Gene,ncol = 4) +
  geom_smooth(method = 'lm',alpha=.3,lwd=.5,col='black',lty=2) +
  annotation_logticks(scaled = TRUE,sides="l") +
  scale_y_log10(breaks = c(0.000001, 10^(-2:5)), 
                labels = c(0, math_format()(-2:5))) +
  scale_x_log10(breaks = c(0.000001, 10^(-2:4)), 
                labels = c(0, math_format()(-2:4))) +
  geom_point(aes(shape = Group,size=3,alpha = .4)) + 
  scale_color_manual(values = colsRS)+
  scale_fill_manual(values = colsRS)+
  theme(legend.position = 'none',text = element_text(size = 14)) + 
  xlab('Faecal Egg Count (eggs/g)') + ylab('Gene-wise transcript count')+
  ggtitle('A')

df2$Gene=substr(df2$Gene,1,13)
dfpl0 = df2[df2$Gene %in% sig,]
dfpl = dfpl0 %>%
  group_by(Gene, Group) %>% 
  summarize(avg_fc = mean(fc), n=n(), sd =sd(fc), se=sd/sqrt(n))

p_qrt = ggplot() +
  geom_point(data = dfpl0,
             aes(x = Group,y = fc), col='black',size = 2, alpha =.5) +
  geom_bar(data = dfpl, aes(x = Group,y = avg_fc,fill = Group,alpha =.6), stat = 'identity') +
  geom_errorbar(data = dfpl, aes(x = Group, ymin = dfpl$avg_fc, ymax = dfpl$avg_fc + dfpl$se,width = 0.05)) +
  facet_wrap(~ Gene, ncol = 4) +
  #scale_color_manual(values = colsRS) +
  scale_fill_manual(values = colsRS) +
  #scale_y_log10() + 
  geom_hline(yintercept = 1) +
  ggtitle('B') +
  ylab('Fold gene relative to S worms') + xlab('Sheep host resistance background') +
  #theme(axis.text.x = element_text(vjust = 0.6, angle = 90),
  theme(text = element_text(size = 14), 
        legend.position = 'none')

pdf(file = '../../../Figure3.pdf',width = 14,height=8)
multiplot(p_rnafec,p_qrt,cols=1)
dev.off()

res2 = NULL
for(i in 1:length(levels(df2$Gene))){
  t = wilcox.test(fc ~ Group,data = df2[df2$Gene==levels(df2$Gene)[i],])$p.value
  tmp2 = data.frame(p=t,g=levels(df2$Gene)[i])
  res2 = rbind(res2,tmp2)
}
res2[res2$p<.05/9,]
#          p             g
# 4  0.00433 HCON_00007272
# 7  0.00433 HCON_00050580
# 11 0.00433 HCON_00182130
# 12 0.00433 HCON_00191700

#### Correlation with RNAseq data
st$Gene = st$gene
mcor = merge(fc.tab[fc.tab$Group == 'R',],st,by='Gene')
mcor = data.frame(mcor)

ggplot(mcor[mcor$Gene!='HCON_00007274',],aes(x = log(fc),y=-log2FoldChange))+
  geom_point(size = 3,alpha =.4) +
  geom_smooth(method='lm') +
  ylab('Fold Change from RNAseq') +
  xlab('Fold Change from qRT-PCR')

rcorr(log(mcor$fc[mcor$Gene!='HCON_00007274']),-mcor$log2FoldChange[mcor$Gene!='HCON_00007274'])
#     x   y
# x 1.0 0.6
# y 0.6 1.0
# 
# n= 5 
# 
# 
# P
#       x     y    
# x       0.286
# y 0.286     

###=============================================
###----------======  WGS =====----------------
###=============================================

require(ggplot2)
theme_set(theme_bw())
require(ggrepel)
require(stringr)
require(reshape2)
require(ggtree)
require(ape)
require(phytools)
## load libraries for the heat map
library("RColorBrewer")
library("gplots")
library("vsn")
require("pheatmap")
library("AnnotationDbi")
library("genefilter")
require("Hmisc")
require("PopGenome")
require("NOISeq") ###-- Diagnostics plot
require(cqn) ###--- Normalization for GC content and gene length if gene x sample bias

setwd('~/Documents/VIR_EVOL/RSG/')

#############--- DNA/RNA extraction 05.10.2017 ---##############

#extract=read.table(file='DNA_RNA_extraction.txt',header=T)
# extract$SheepID=factor(extract$SheepID)
# #- DNA
# ggplot(extract,aes(x=SheepID,y=DNAyield,color=Group))+
#   geom_point() +
#   theme_bw()
# a=aggregate(DNAyield ~ Group, FUN=mean,data=extract)
# a$std=aggregate(DNAyield ~ Group, FUN=sd,data=extract)
# # Group DNAyield std.Group std.DNAyield
# # 1     R 390.5263         R     192.6457
# # 2     S 678.8000         S     294.7751
# kruskal.test(DNAyield ~ Group,data=extract)
# # 
# # Kruskal-Wallis rank sum test
# # 
# # data:  DNAyield by Group
# # Kruskal-Wallis chi-squared = 10.262, df = 1, p-value = 0.001358
# 
# ggplot(extract,aes(x=SheepID,y=DNAyield,fill=Group))+
#   geom_boxplot() +
#   theme_bw()
# ggplot(extract,aes(x=Group,y=DNAyield,color=Group))+
#   geom_point() +
#   theme_bw()
# 
# #- RNA
# ggplot(extract,aes(x=SheepID,y=RNAyield,color=Group))+
#   geom_point() +
#   scale_y_log10() +
#   theme_bw()
# a=aggregate(RNAyield ~ Group, FUN=mean,data=extract)
# a$std=aggregate(RNAyield ~ Group, FUN=sd,data=extract)
# #   Group  RNAyield std.Group std.RNAyield
# # 1     R  511.2274         R     717.2874
# # 2     S 1033.8740         S    1687.3028
# kruskal.test(RNAyield ~ Group,data=extract)
# # Kruskal-Wallis rank sum test
# # 
# # data:  RNAyield by Group
# # Kruskal-Wallis chi-squared = 0.087057, df = 1, p-value = 0.768
# 
# ggplot(extract,aes(x=SheepID,y=RNAyield,fill=Group))+
#   geom_boxplot() +
#   theme_bw()
# 
# rnaQ=extract[extract$RNAyield>100,]
# table(rnaQ$Group)
# # R  S 
# # 12 14 
# 
# ggplot(rnaQ,aes(x=SheepID,y=RNAyield,color=Group))+
#   geom_point() +
#   theme_bw()
# 
# qq=extract[extract$RNAyield>=50,]
# ggplot(extract,aes(x=DNAyield,y=RNAyield,color=Group))+
#   geom_point()+
#   theme_bw()

#############--- ANALYSIS DNA ---##############

###--- FST
setwd('/Users/gsalle/Documents/INRA/GEMANEMA/Paper/data/WGS/')

chromlist = c(1,2,3,4,5,'X')
prefix = 'hcontortus_chr'
suffix = '_Celeg_TT_arrow_pilon'

##-- 1R 2S 3R 4S 5S 6R 7S 8R
#-- Min coverage has min impact on the results
fst1 = read.table(file='./rsg.tot.fst',
                  header=F,sep='\t')

#-- vectors to indicate whether Fst contrast is within 1 or betw 2 exp groups
#- 1 = R-R / 2 = R-S / 3 = S-S
comp=c(2,1,2,2,1,2,1,
       2,3,3,2,3,2,
       2,2,1,2,1,
       3,2,3,2,
       2,3,2,
       2,1,
       2)

#-- convert
for(i in 6:(length(comp)+5)){
  temp = as.character(fst1[,i])
  fst1[,i] = as.numeric(str_split_fixed(fst1[,i],"=",2)[,2])
}

#-- Coverage found
summary(fst1$V5)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 20.00   39.80   45.70   46.08   51.80  209.40 

quantile(fst1$V5)
#   0%   25%   50%   75%  100% 
# 20.0  39.8  45.7  51.8 209.4 

quantile(fst1$V4)
#    0%   25%   50%   75%  100% 
# 0.011 0.505 0.772 0.893 1.000

summary(fst1$V3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.0   248.0   411.0   423.1   574.0  1638.0 

## Remove windows with less than 50 SNPs or with too high coverage
fst1 = fst1[fst1$V3 > 100 & fst1$V5 < 120 & fst1$V4 > 0.5,]

#####-------------------------------------------- split by within/between
chrom.colors = rev((viridis_pal(option='D'))(6))
#chrom.colors=c(cols[1],cols[12],cols[2],cols[3],cols[20],cols[21])

fst1$fstb = rowMeans(fst1[,5+which(comp==2)])
fst1$fstwR = rowMeans((fst1[,5+which(comp==1)]))
fst1$fstwS = rowMeans(fst1[,5+which(comp==3)])
fst1$fstw = rowMeans(fst1[,5+which(comp==3 |comp==1)])
fst1$bin = seq(1:dim(fst1)[1])
fst = melt(fst1[,c('V1','V2','V3','V4','V5','bin','fstb','fstwR','fstwS','fstw')],1:6)


#### Suppl Fig 6
fstplot = fst[fst$variable %in% c('fstb','fstw'),]
fstplot$variable = as.character(fstplot$variable)
fstplot$variable[fstplot$variable=='fstb']='FST across sheep genotype'
fstplot$variable[fstplot$variable=='fstw']='FST within sheep genotype'
fstplot$variable=factor(fstplot$variable)

# retrieve genomic position
chrom_end = c(aggregate(V2 ~ V1, FUN=max,data=fstplot)$V2)
chrom_end = c(0,chrom_end[-6])
last=0
gp = NULL
for(i in 1:length(chrom_end)){
  gp[i] = sum(chrom_end[1:i])
}
gp

fstplot$genpos = fstplot$V2 + gp[match(fstplot$V1,levels(fstplot$V1))]

## Genome-wide plot : no diff between /within
tiff(file = paste0(draft,'supplementary_Figure6.tif'),units='in', width = 14, height = 8,res=600)

png(file = paste0(draft,'supplementary_Figure6.png'),
    width = 14, height = 8,units='in',res =600)

ggplot(fstplot,aes(x=genpos/1e6,y=value,col=V1)) +
  scale_y_continuous(limits = c(0,.5),breaks = seq(0,.5,0.1)) +
  scale_x_continuous(limits = c(0,285),breaks = seq(0,285,15)) +
  scale_color_manual(values=c(chrom.colors)) +
  geom_point(size=.1,alpha = .3) +
  facet_wrap(~ variable,ncol=1) + theme_bw() +
  #geom_hline(yintercept = mean(fstplot$value) + 3*sd(fstplot$value)) +
  xlab('Genomic Position (Mbp)') +
  ylab('Mean pair-wise FST') +
  theme(legend.position = 'none',text = element_text(size=16))
dev.off()

## Summary
summary(fst1$fstb)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.003708 0.013831 0.018214 0.018567 0.022531 0.080006 

summary(fst1$fstw)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.003708 0.013907 0.018345 0.018705 0.022717 0.090773 

summary(fst1$fstwR)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.003471 0.013431 0.017817 0.018551 0.022536 0.097233

summary(fst1$fstwS)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.002814 0.013532 0.017957 0.018860 0.022825 0.091378 

## Another representation
ggplot(fst1,aes(x = bin, y=fstb-fstw,col=V1)) + 
  geom_point(size=.1) + 
  scale_color_manual(values = chrom.colors) +
  theme(legend.position = 'none',text=element_text(size=16))

mean(fst1$fstb-fst1$fstw)
#[1] -0.0001372572

## FST distribution showing range of values and lack of difference
require(ggridges)
levels(fst$variable) = list(FST_between='fstb',FST_within_R='fstwR',FST_within_S='fstwS',FST_within='fstw')

#pdf(file = paste0(draft,'supplementary_Figure7.pdf'), width = 14, height = 8)
# ggplot(fst,aes(x = value,y = factor(variable),fill = factor(variable))) + 
#   geom_density_ridges(scale = 0.9) + 
#   theme(legend.position = 'bottom', text = element_text(size=16)) +
#   ylab('') + xlab('FST estimates') + theme(legend.position = 'none')
#dev.off()

### Plot within R / within S
# pdf(file = paste0(draft,'supplementary_Figure6.pdf'), width = 14, height = 8)
# ggplot(fst,aes(x=genpos/1e6,y=value,col=V1)) +
#   scale_y_continuous(limits = c(0,.5),breaks = seq(0,.5,0.1)) +
#   scale_x_continuous(limits = c(0,285),breaks = seq(0,285,15)) +
#   scale_color_manual(values=c(chrom.colors)) +
#   geom_point(size=.1,alpha = .3) +
#   facet_wrap(~ variable,ncol=1) + theme_bw() +
#   #geom_hline(yintercept = mean(fstplot$value) + 3*sd(fstplot$value)) +
#   xlab('Genomic Position (Mbp)') +
#   ylab('Mean pair-wise FST') +
#   theme(legend.position = 'none',text = element_text(size=16))
# dev.off()

## Outlier regions between groups that are not outlier in within-group comparison
cutoffb = mean(fst$value[fst$variable=='FST_between']) + 3*sd(fst$value[fst$variable=='FST_between'])
cutoffw = mean(fst$value[fst$variable=='FST_within']) + 3*sd(fst$value[fst$variable=='FST_within'])

outliers = fst1[fst1$fstb > cutoffb & fst1$fstw < cutoffw,]
table(outliers$V1)
# hcontortus_chr1_Celeg_TT_arrow_pilon hcontortus_chr2_Celeg_TT_arrow_pilon 
#                                   34                                   82 
# hcontortus_chr3_Celeg_TT_arrow_pilon hcontortus_chr4_Celeg_TT_arrow_pilon 
#                                   30                                   67 
# hcontortus_chr5_Celeg_TT_arrow_pilon hcontortus_chrX_Celeg_TT_arrow_pilon 
#                                   95                                  171 
fstplot$outliers = 0
fstplot$outliers[which(paste0(fstplot$V1,fstplot$V2) %in% paste0(outliers$V1,outliers$V2))] = 1

aggregate(value ~ outliers,data=fstplot,FUN=mean)
#   outliers      value
# 1        0 0.01682341
# 2        1 0.03269686

# load the gtf file
gtf.gr <- rtracklayer::import(con ="~/Documents/REF_FILES/haemonchus_contortus.PRJEB506.WBPS13.annotations.gtf.gz", format = "gtf" )
gtf.gr <- sortSeqlevels(gtf.gr) # make sure it is sorted
gtf.gr <- sort(gtf.gr)

## Convert file into bed file
out = outliers[,c('V1','V2')]
bed.candid2 = with(out, GRanges(V1, IRanges(start = V2, end = V2)))

hits = findOverlaps(gtf.gr, bed.candid2)
length(unique(gtf.gr$gene_id[queryHits(hits)]))
#[1] 37
factor(substr(unique(gtf.gr$gene_id[queryHits(hits)]),6,18))
# [1] HCON_00022360 HCON_00025640 HCON_00036120 HCON_00049130 HCON_00049880 HCON_00193140
# [7] HCON_00052990 HCON_00094050 HCON_00091320 HCON_00094060 HCON_00121240 HCON_00128300
# [13] HCON_00109020 HCON_00141740 HCON_00143410 HCON_00143420 HCON_00136900 HCON_00140940
# [19] HCON_00146000 HCON_00164360 HCON_00168080 HCON_00168280 HCON_00169420 HCON_00185600
# [25] HCON_00186210 HCON_00190240 HCON_00190250 HCON_00165005 HCON_00165250 HCON_00166090
# [31] HCON_00168010 HCON_00169480 HCON_00184620 HCON_00185700 HCON_00186770 HCON_00188290
# [37] HCON_00188540

## Any DE genes in these ?
dfin=c('HCON_00007272','HCON_00007274',
       'HCON_00029050','HCON_00050580',
       'HCON_00087240','HCON_00090580',
       'HCON_00132830','HCON_00182130','HCON_00191700')
match(dfin,factor(sapply(str_split(unique(gtf.gr$gene_id[queryHits(hits)]),'gene:'),function(x) x[2])))
#[1] NA NA NA NA NA NA NA NA NA NA

### Show that outlier regions are identical
fstpl2 = fstplot[fstplot$outliers==1,]
fstpl2$value[fstpl2$variable=='FST.within.Sheep.Genotype'] = -1*fstpl2$value[fstpl2$variable=='FST.within.Sheep.Genotype']

ggplot(fstpl2,aes(x=bin,y=value,col=V1)) +
  scale_color_manual(values=c(chrom.colors)) +
  geom_point(size=.8) +
  #facet_wrap(~ variable,ncol=1) + theme_bw() +
  #geom_hline(yintercept = mean(fstplot$value) + 3*sd(fstplot$value)) +
  xlab('bin')+
  ylab('Mean pair-wise FST') +
  theme(legend.position = 'none',text=element_text(size=16))



## Find closest match between comparisons
ofstb = fst1[fst1$fstb > cutoffb,]
ofstw = fst1[fst1$fstw > cutoffw,]
dim(ofstb)
#[1] 1717   38
dim(ofstw)
#[1] 2018   38
distbtow = array(-1,dim(ofstb)[1])
distwtob = array(-1,dim(ofstw)[1])

for(i in 1:dim(ofstb)[1]){
  chr = ofstb$V1[i]
  pos = ofstb$V2[i]
  distbtow[i] = abs(ofstw$V2[ofstw$V1==chr]-pos)[which.min(abs(ofstw$V2[ofstw$V1==chr]-pos))]
}
for(i in 1:dim(ofstw)[1]){
  chr = ofstw$V1[i]
  pos = ofstw$V2[i]
  distwtob[i] = abs(ofstb$V2[ofstb$V1==chr]-pos)[which.min(abs(ofstb$V2[ofstb$V1==chr]-pos))]
}
plot(density(distwtob))
plot(density(distbtow))

summary(distwtob)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0       0       0  118501    7000 9338000
summary(distbtow)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0       0       0  122531   26000 4389000

####------------ FST-based tree between samples
library(ape)
vecmfst=as.numeric(colMeans(fst1[,6:dim(fst1)[2]]))
Fstoverall<-matrix(ncol=8,nrow=8)
n=0
for (j in 1:7) {
  for (i in (j+1):8) {
    n=n+1
    Fstoverall[i,j]=vecmfst[n]
  }
}
sn = c('R_1','S_1','R_2','S_2','S_3','R_3','S_4','R_4')
par(mfrow=c(1,1))
Fstoverall2 = Fstoverall/(1-Fstoverall)
inbreed = matrix(Fstoverall2,ncol=8,dimnames=list(sn[1:8],sn[1:8]))
dd2.1 = as.dist(inbreed)
require(ggtree)
a = nj(dd2.1)
groupInfo <- split(a$tip.label, gsub("_\\w+", "", a$tip.label))
a <- groupOTU(a, groupInfo)

ggtree(a, aes(color=group), layout = 'rectangular') + 
  geom_tiplab(cex=5) +
  geom_treescale(fontsize = 4,linesize = 1) +
  theme(title=element_text(size=18))+
  ggtitle("NJ tree based on pair-wise FST")

sessionInfo()
# R version 3.5.0 (2018-04-23)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.14.5
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] splines   grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] scales_1.0.0                DESeq2_1.22.2               edgeR_3.24.3                limma_3.38.3               
# [5] GenomicFeatures_1.34.8      tximportData_1.10.0         readr_1.3.1                 tximport_1.10.1            
# [9] rhdf5_2.26.2                sleuth_0.30.0               cqn_1.28.1                  quantreg_5.38              
# [13] preprocessCore_1.44.0       nor1mix_1.2-3               mclust_5.4.3                NOISeq_2.26.1              
# [17] Matrix_1.2-17               PopGenome_2.7.1             ff_2.2-14                   bit_1.1-14                 
# [21] Hmisc_4.2-0                 Formula_1.2-3               survival_2.44-1.1           lattice_0.20-38            
# [25] VennDiagram_1.6.20          futile.logger_1.4.3         SummarizedExperiment_1.12.0 DelayedArray_0.8.0         
# [29] BiocParallel_1.16.6         matrixStats_0.54.0          GenomicRanges_1.34.0        GenomeInfoDb_1.18.2        
# [33] topGO_2.34.0                SparseM_1.77                GO.db_3.7.0                 graph_1.60.0               
# [37] genefilter_1.64.0           AnnotationDbi_1.44.0        IRanges_2.16.0              S4Vectors_0.20.1           
# [41] pheatmap_1.0.12             Biobase_2.42.0              BiocGenerics_0.28.0         gplots_3.0.1.1             
# [45] RColorBrewer_1.1-2          phytools_0.6-60             maps_3.3.0                  ape_5.3                    
# [49] ggtree_1.14.6               reshape2_1.4.3              stringr_1.4.0               ggrepel_0.8.1              
# [53] ggplot2_3.1.1              
# 
# loaded via a namespace (and not attached):
# [1] backports_1.1.4          fastmatch_1.1-0          plyr_1.8.4               igraph_1.2.4.1           lazyeval_0.2.2          
# [6] digest_0.6.19            htmltools_0.3.6          gdata_2.18.0             magrittr_1.5             checkmate_1.9.3         
# [11] memoise_1.1.0            cluster_2.0.9            Biostrings_2.50.2        annotate_1.60.1          prettyunits_1.0.2       
# [16] colorspace_1.4-1         blob_1.1.1               xfun_0.7                 dplyr_0.8.1              crayon_1.3.4            
# [21] RCurl_1.95-4.12          jsonlite_1.6             phangorn_2.5.3           glue_1.3.1               gtable_0.3.0            
# [26] zlibbioc_1.28.0          XVector_0.22.0           MatrixModels_0.4-1       Rhdf5lib_1.4.3           futile.options_1.0.1    
# [31] DBI_1.0.0                Rcpp_1.0.1               plotrix_3.7-5            progress_1.2.2           xtable_1.8-4            
# [36] htmlTable_1.13.1         tidytree_0.2.4           foreign_0.8-71           animation_2.6            httr_1.4.0              
# [41] htmlwidgets_1.3          acepack_1.4.1            pkgconfig_2.0.2          XML_3.98-1.19            nnet_7.3-12             
# [46] locfit_1.5-9.1           labeling_0.3             tidyselect_0.2.5         rlang_0.3.4              munsell_0.5.0           
# [51] tools_3.5.0              RSQLite_2.1.1            yaml_2.2.0               knitr_1.23               bit64_0.9-7             
# [56] caTools_1.17.1.2         purrr_0.3.2              nlme_3.1-140             formatR_1.6              biomaRt_2.38.0          
# [61] compiler_3.5.0           rstudioapi_0.10          affyio_1.52.0            treeio_1.6.2             clusterGeneration_1.3.4 
# [66] tibble_2.1.1             geneplotter_1.60.0       stringi_1.4.3            pillar_1.4.1             BiocManager_1.30.4      
# [71] combinat_0.0-8           data.table_1.12.2        bitops_1.0-6             rtracklayer_1.42.2       R6_2.4.0                
# [76] latticeExtra_0.6-28      affy_1.60.0              KernSmooth_2.23-15       gridExtra_2.3            lambda.r_1.2.3          
# [81] MASS_7.3-51.4            gtools_3.8.1             assertthat_0.2.1         withr_2.1.2              GenomicAlignments_1.18.1
# [86] Rsamtools_1.34.1         mnormt_1.5-5             GenomeInfoDbData_1.2.0   expm_0.999-4             hms_0.4.2               
# [91] quadprog_1.5-7           rpart_4.1-15             tidyr_0.8.3              coda_0.19-2              rvcheck_0.1.3           
# [96] numDeriv_2016.8-1        scatterplot3d_0.3-41     base64enc_0.1-3 