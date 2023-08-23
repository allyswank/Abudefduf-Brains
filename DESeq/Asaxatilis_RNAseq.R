#############################################################################
########################## RNA-seq Gene Expression ##########################
###################### Abudefduf saxatilis - Full Brain #####################
#############################################################################

#Load required packages
library(ggplot2)
library(DESeq2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(tidyr)

#Set working directory and load the raw gene counts table and metadata file
#setwd("C:/Users/allys/Box/Auburn/Bernal Lab/Thesis/Data")
setwd("C:/Users/allys/Box/Auburn/Scientific Writing/Analyses and Data/DESeq")

gene_counts <- as.matrix(read.table("allcounts.txt", row.names = 1, header=TRUE, sep = "\t"))
gene_counts <- gene_counts[,colnames(gene_counts)!="X.1"]

metadata = read.csv("a.saxatilis_attributes.csv", header = TRUE, sep = ",", row.names = 1)
metadata$Temperature = as.factor(metadata$Temperature)
metadata$Complexity = as.factor(metadata$Complexity)
metadata$Placement = as.factor(metadata$Placement)
metadata$Tower = as.factor(metadata$Tower)
metadata$Treatment = as.factor(metadata$Treatment)

#making sure sample ID's are the same in raw counts and metadata files
#This should be "TRUE"
all(rownames(metadata) == colnames(gene_counts))


#############################################################################
####################### Sample x Sample Differences? ########################
#############################################################################
dds <- DESeqDataSetFromMatrix(countData = gene_counts, colData = metadata, design =~ Treatment)  
dds <- DESeq(dds)

vsd <- getVarianceStabilizedData(dds)
vsd.values <- cbind(results$pvalue, results$padj)
colnames(vsd.values) <- c("pval", "padj")
vsd.pvalues <- cbind(vsd, vsd.values)

vsdata <- vst(dds, blind=FALSE)
sampleDists <- dist(t(vsd))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsdata))
colnames(sampleDistMatrix) <- paste(colnames(vsdata))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(250)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows=sampleDists, 
         clustering_distance_cols=sampleDists, 
         col=colors)
#There are two outliers here that could be skewing results, so they will be removed and analyses will be rerun

############################################################################

gene_counts <- gene_counts[,colnames(gene_counts) != "C.HC.09"]
gene_counts <- gene_counts[,colnames(gene_counts) != "HW.HC.02"]
metadata <- metadata[rownames(metadata) != "C.HC.09",]
metadata <- metadata[rownames(metadata) != "HW.HC.02",]

all(rownames(metadata) == colnames(gene_counts))

colMeans(gene_counts)

#############################################################################

#PCA of all fish based on sig treatment genes

dds <- DESeqDataSetFromMatrix(countData = gene_counts, colData = metadata, design =~ Treatment + Placement)  
dds <- DESeq(dds)  
res <- results(dds, alpha = 0.05)
resultsNames(dds)

LRT_All <- DESeq(dds, test = "LRT", reduced =~ Placement)
res_LRT_All <- results(LRT_All, alpha = 0.05)
summary(res_LRT_All)
#out of 117945 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 171, 0.14%
#LFC < 0 (down)     : 114, 0.097%
#outliers [1]       : 8, 0.0068%
#low counts [2]     : 102807, 87%

table(res_LRT_All$padj<0.05)
#FALSE  TRUE 
#14845   285  genes that are different based on the specific individual conditions 
vsdata <- vst(dds)
sigs_all <- subset(res_LRT_All, padj < 0.05) #We want just the significant genes, so isolate them here
summary(sigs_all)
sigs_names_all <- rownames(sigs_all) #Isolate just the names of the significant genes
VST.all.subset <- vsdata[rownames(vsdata) %in% sigs_names_all, ] #Extract transformed counts for significant genes
summary(VST.all.subset) #Make sure you've got them and they are correct
plotPCA(VST.all.subset, intgroup="Treatment")
pca <- plotPCA(VST.all.subset, intgroup="Treatment", ntop = nrow(VST.all.subset), returnData = TRUE) #effect of individual traits

write.csv(res_LRT_All, file="LRT_All.csv", quote=FALSE)


pca_plot = ggplot(data=pca, aes_string(x="PC1", y="PC2", color="Treatment")) + 
  geom_point(size=4, alpha = 0.6, stroke = NA) + 
  xlab(paste0("PC1: 36% variance")) +
  ylab(paste0("PC2: 14% variance")) +
  coord_fixed() +
  scale_color_manual(values=c("#009E73", "#56B4E9", "#E69F00", "#D55E00")) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 15, colour = "black"))
pca_plot

pdf("PCA_AllSamples.pdf")
pca_plot
dev.off()

#Likelihood Ratio Tests (LRT)
#comparing the effects of one variable on all groups together

dds <- DESeqDataSetFromMatrix(countData = gene_counts, colData = metadata, design =~ Temperature + Complexity)
dds <- DESeq(dds)
results <- results(dds)
resultsNames(dds)


#LowComplexity vs. HighComplexity
LRT_LC_HC <- DESeq(dds, test = "LRT", reduced =~ Temperature)
res_LRT_LC_HC <- results(LRT_LC_HC, alpha = 0.05)
summary(res_LRT_LC_HC)
  #No DEGs between LC and HC (including both Temperature groups)


#Heatwave vs. Control
LRT_C_HW <- DESeq(dds, test = "LRT", reduced =~ Complexity)
res_LRT_C_HW <- results(LRT_C_HW, alpha = 0.05)
summary(res_LRT_C_HW)

  #Total of 255 DEGs between HW and C temperatures (including both Complexity groups)
C_HW_LRTvst <- vst(LRT_C_HW) #Transform the count data using variance stabilizing transformation
sigs_LRT_C_HW <- subset(res_LRT_C_HW, padj < 0.05) #We want just the significant genes, so isolate them here
sigs_names_C_HW <- rownames(sigs_LRT_C_HW) #Isolate just the names of the significant genes
VST.C_HW.subset <- C_HW_LRTvst[rownames(C_HW_LRTvst) %in% sigs_names_C_HW, ] #Extract transformed counts for significant genes
summary(VST.C_HW.subset) #Make sure you've got them and they are correct


#Plot PCA of significant genes associated with temperature treatments (includes both Complexity groups)
plotPCA(VST.C_HW.subset, intgroup = "Temperature")

write.csv(res_LRT_C_HW, file="LRT_CvsHW.csv", quote=FALSE)
#measures <- read.csv('LRT_CvsHW.csv', header = TRUE, sep = ",")
#measures <- na.omit(measures)
#write.csv(measures, file="LRT_CvsHW.csv", quote=FALSE)

#############################################################################

#Wald test for pairwise comparisons

#Control temperatures: HighComplexity vs LowComplexity
metadata_C_HCvsLC <- subset(metadata, Temperature == "27")
genecounts_C_HCvsLC <- gene_counts[,c(1:19)]
metadata_C_HCvsLC$Complexity = as.factor(metadata_C_HCvsLC$Complexity)

Cdds <- DESeqDataSetFromMatrix(countData = genecounts_C_HCvsLC, colData = metadata_C_HCvsLC, design =~ Complexity)
Cdds <- DESeq(Cdds)
resultsNames(Cdds)

C_HCvsLC <- results(Cdds, alpha = 0.05, contrast = c("Complexity", "L", "H"))
C_HCvsLC
sigs_C_HCvsLC <- subset(C_HCvsLC, padj < 0.05)
summary(C_HCvsLC)
  #No significant DEGs for Complexity in Control temp groups



#Heatwave temperatures: HighComplexity vs LowComplexity
metadata_HW_HCvsLC <- subset(metadata, Temperature == "31")
genecounts_HW_HCvsLC <- gene_counts[,c(20:38)]
metadata_HW_HCvsLC$Complexity = as.factor(metadata_HW_HCvsLC$Complexity)

HWdds <- DESeqDataSetFromMatrix(countData = genecounts_HW_HCvsLC, colData = metadata_HW_HCvsLC, design =~ Complexity)
HWdds <- DESeq(HWdds)
resultsNames(HWdds)

HW_HCvsLC <- results(HWdds, alpha = 0.05, contrast = c("Complexity", "L", "H"))
HW_HCvsLC
sigs_HW_HCvsLC <- subset(HW_HCvsLC, padj < 0.05)
summary(HW_HCvsLC)
  #No significant DEGs for Complexity in Heatwave temp groups



#HighComplexity: Control vs. Heatwave Temps
metadata_HC_CvsHW <- subset(metadata, Complexity == "H")
genecounts_HC_CvsHW <- gene_counts[,c(1:9,20:28)]
metadata_HC_CvsHW$Temperature = as.factor(metadata_HC_CvsHW$Temperature)

HCdds <- DESeqDataSetFromMatrix(countData = genecounts_HC_CvsHW, colData = metadata_HC_CvsHW, design =~ Temperature)
HCdds <- DESeq(HCdds)
resultsNames(HCdds)

HC_CvsHW <- results(HCdds, alpha = 0.05, contrast = c("Temperature", "31", "27"))
HC_CvsHW
sigs_HC_CvsHW <- subset(HC_CvsHW, padj < 0.05)
summary(HC_CvsHW)
  #Total of 56 DEGs between Control and Heatwave Temps when habitat complexity was high

vst_HC_CvsHW <- vst(HCdds)
signames_HC_CvsHW <- rownames(sigs_HC_CvsHW)
vst.HC_CvsHW.subset <- vst_HC_CvsHW[rownames(vst_HC_CvsHW) %in% signames_HC_CvsHW, ]
summary(vst.HC_CvsHW.subset)
plotPCA(vst.HC_CvsHW.subset, intgroup = "Temperature")

write.csv(HC_CvsHW, file="Wald_HC_CvsHW.csv", quote=FALSE)
#measures <- read.csv('Wald_HC_CvsHW.csv', header = TRUE, sep = ",")
#measures <- na.omit(measures)
#write.csv(measures, file="Wald_HC_CvsHW.csv", quote=FALSE)

#LowComplexity: Control vs. Heatwave Temps
metadata_LC_CvsHW <- subset(metadata, Complexity == "L")
genecounts_LC_CvsHW <- gene_counts[,c(10:19,29:38)]
metadata_LC_CvsHW$Temperature = as.factor(metadata_LC_CvsHW$Temperature)

LCdds <- DESeqDataSetFromMatrix(countData = genecounts_LC_CvsHW, colData = metadata_LC_CvsHW, design =~ Temperature)
LCdds <- DESeq(LCdds)
resultsNames(LCdds)

LC_CvsHW <- results(LCdds, alpha = 0.05, contrast = c("Temperature", "31", "27"))
LC_CvsHW
sigs_LC_CvsHW <- subset(LC_CvsHW, padj < 0.05)
summary(LC_CvsHW)
  #Total of 122 DEGs between Control and Heatwave Temps when habitat complexity was low

vst_LC_CvsHW <- vst(LCdds)
signames_LC_CvsHW <- rownames(sigs_LC_CvsHW)
vst.LC_CvsHW.subset <- vst_LC_CvsHW[rownames(vst_LC_CvsHW) %in% signames_LC_CvsHW, ]
summary(vst.LC_CvsHW.subset)
plotPCA(vst.LC_CvsHW.subset, intgroup = "Temperature")

write.csv(LC_CvsHW, file="Wald_LC_CvsHW.csv", quote=FALSE)
#measures <- read.csv('Wald_LC_CvsHW.csv', header = TRUE, sep = ",")
#measures <- na.omit(measures)
#write.csv(measures, file="Wald_LC_CvsHW.csv", quote=FALSE)

#LowComplexity Control vs. HighComplexity Heatwave
metadata_CLCvsHWHC <- metadata[c(10:28),]
genecounts_CLCvsHWHC <- gene_counts[,c(10:28)]
metadata_CLCvsHWHC$Treatment = as.factor(metadata_CLCvsHWHC$Treatment)

dds <- DESeqDataSetFromMatrix(countData = genecounts_CLCvsHWHC, colData = metadata_CLCvsHWHC, design =~ Treatment)
dds <- DESeq(dds)
resultsNames(dds)

CLCvsHWHC <- results(dds, alpha = 0.05, contrast = c("Treatment", "HWHC", "CLC"))
CLCvsHWHC
sigs_CLCvsHWHC <- subset(CLCvsHWHC, padj < 0.05)
summary(CLCvsHWHC)
  #Total of 102 DEGs

vst_CLCvsHWHC <- vst(dds)
signames_CLCvsHWHC <- rownames(sigs_CLCvsHWHC)
vst.CLCvsHWHC.subset <- vst_CLCvsHWHC[rownames(vst_CLCvsHWHC) %in% signames_CLCvsHWHC, ]
summary(vst.CLCvsHWHC.subset)
plotPCA(vst.CLCvsHWHC.subset, intgroup = "Treatment")

write.csv(CLCvsHWHC, file="Wald_CLCvsHWHC.csv", quote=FALSE)
#measures <- read.csv('Wald_CLCvsHWHC.csv', header = TRUE, sep = ",")
#measures <- na.omit(measures)
#write.csv(measures, file="Wald_CLCvsHWHC.csv", quote=FALSE)


#HighComplexity Control vs. LowComplexity Heatwave
metadata_CHCvsHWLC <- metadata[c(1:9,29:38),]
genecounts_CHCvsHWLC <- gene_counts[,c(1:9,29:38)]
metadata_CHCvsHWLC$Treatment = as.factor(metadata_CHCvsHWLC$Treatment)

dds <- DESeqDataSetFromMatrix(countData = genecounts_CHCvsHWLC, colData = metadata_CHCvsHWLC, design =~ Treatment)
dds <- DESeq(dds)
resultsNames(dds)

CHCvsHWLC <- results(dds, alpha = 0.05, contrast = c("Treatment", "HWLC", "CHC"))
CHCvsHWLC
sigs_CHCvsHWLC <- subset(CHCvsHWLC, padj < 0.05)
summary(CHCvsHWLC)
  #Total of 96 DEGs

vst_CHCvsHWLC <- vst(dds)
signames_CHCvsHWLC <- rownames(sigs_CHCvsHWLC)
vst.CHCvsHWLC.subset <- vst_CLCvsHWHC[rownames(vst_CHCvsHWLC) %in% signames_CHCvsHWLC, ]
summary(vst.CHCvsHWLC.subset)
plotPCA(vst.CHCvsHWLC.subset, intgroup = "Treatment")

write.csv(CHCvsHWLC, file="Wald_CHCvsHWLC.csv", quote=FALSE)
#measures <- read.csv('Wald_CHCvsHWLC.csv', header = TRUE, sep = ",")
#measures <- na.omit(measures)
#write.csv(measures, file="Wald_CHCvsHWLC.csv", quote=FALSE)

#############################################################################
############################## Random Effects? ###############################
#############################################################################

Tdds <- DESeqDataSetFromMatrix(countData = gene_counts, colData = metadata, design =~ Tower + Placement)
Tdds <- DESeq(Tdds)
resultsNames(Tdds)


Tower <- DESeq(Tdds, test = "LRT", reduced =~ Placement)
res_Tower <- results(Tower, alpha = 0.05)
summary(res_Tower)
  #Total of 424 DEGs between all towers
Tower_LRTvst <- vst(Tower) #Transform the count data using variance stabilizing transformation
sigs_Tower <- subset(res_Tower, padj < 0.05) #We want just the significant genes, so isolate them here
sigs_names_Tower <- rownames(sigs_Tower) #Isolate just the names of the significant genes
VST.Tower.subset <- Tower_LRTvst[rownames(Tower_LRTvst) %in% sigs_names_Tower, ] #Extract transformed counts for significant genes
summary(VST.Tower.subset) #Make sure you've got them and they are correct

plotPCA(VST.Tower.subset, intgroup = "Tower")


Position <- DESeq(Tdds, test = "LRT", reduced =~ Tower)
res_Position <- results(Position, alpha = 0.05)
summary(res_Position)
  #No DEGs between position in the tower



############################################################################
################################# FIGURES ##################################
############################################################################

#install.packages("VennDiagram")
library(VennDiagram)

venn.diagram(x = list(signames_HC_CvsHW, signames_LC_CvsHW), 
             category.names = c("HC_CvsHW", "LC_CvsHW"), 
             filename = 'VennDiagram-HC_CvsHW-LC_CvsHW.png', output = TRUE)

venn.diagram(x = list(signames_HC_CvsHW, signames_LC_CvsHW, signames_CHCvsHWLC, signames_CLCvsHWHC), 
             category.names = c("HC_CvsHW", "LC_CvsHW", "CHCvsHWLC", "CLCvsHWHC"), 
             filename = 'VennDiagram-all-sig-groups.png', 
             output = TRUE)

#Finding the shared genes between comparisons! 

#HighComplexity CvsHW - LowComplexity CvsHW

tempgenes <- intersect(sigs_names_C_HW, sigs_names_C_HW)
tempgenes
write.csv(tempgenes, file="TempGenesLRT.csv", quote=FALSE)


shared.HC_CvsHW.LC_CvsHW <- intersect(signames_HC_CvsHW, signames_LC_CvsHW)
shared.HC_CvsHW.LC_CvsHW
write.csv(shared.HC_CvsHW.LC_CvsHW, file="HeatwaveResponseGenes.csv", quote=FALSE)

LC <- as.data.frame(sigs_LC_CvsHW)
HC <- as.data.frame(sigs_HC_CvsHW)

#isolating the private genes associated with complexity
all <- union(signames_HC_CvsHW, signames_LC_CvsHW)
non.matched <- all[!all %in% shared.HC_CvsHW.LC_CvsHW]
LC.genes <- intersect(non.matched, signames_LC_CvsHW)
LC.genes
HC.genes <- intersect(non.matched, signames_HC_CvsHW)
HC.genes
write.csv(LC.genes, file="LC_genes.csv", quote=FALSE)
write.csv(HC.genes, file = "HC_genes.csv", quote = FALSE)

LCpriv <- as.data.frame(LC.genes, row.names = LC.genes)
#write.csv(LCpriv, file = "LC_genes.csv", quote = FALSE)
HCpriv <- as.data.frame(HC.genes, row.names = HC.genes)

LC_private <- merge(LC, LCpriv, by = 'row.names')
HC_private <- merge(HC, HCpriv, by = 'row.names')


shared.all.sigs <- intersect(signames_CHCvsHWLC, signames_CLCvsHWHC) #33 genes shared here
write.csv(shared.all.sigs, file = "shared_CHC-HWLC&CLC-HWHC.csv", quote = FALSE)
shared.all.sigs <- intersect(shared.all.sigs, signames_HC_CvsHW)
shared.all.sigs <- intersect(shared.all.sigs, signames_LC_CvsHW)
shared.all.sigs

