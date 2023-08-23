#############################################################################
############################## WGCNA - RNA-seq ##############################
##################### Abudefduf saxatilis - Full Brains #####################
#############################################################################

#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-RelateModsToTraits.pdf

#this is the WGCNA for Abudefduf saxatilis - Thermal and Habitat stress 

#installing WGCNA again
#install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
#install.packages("BiocManager")
#BiocManager::install("WGCNA")
#BiocManager::install("DESeq2")

setwd("C:/Users/allys/Box/Bernal_lab/Ally/Chapter1/Analyses/WGCNA")
library('DESeq2')
library('pheatmap')
library('RColorBrewer')
library('ggplot2')
library('WGCNA')
library('tidyverse')

options(stringsAsFactors = FALSE)

x<-read.csv('allcounts.txt', row.names="Gene", header=T, sep='\t')
x <- x[,colnames(x)!="X"]

head(x)
nrow(x)#119863


#filter low counts
x$filtr<-apply(x, 1, function(k) mean(k > 10)) > 0.5
x<-x[!(x$filtr=="FALSE"),]
nrow(x) #10658

#This is the code to only analyze the 50% with the highest variation
#counts$variance = apply(counts, 1, var)
#counts2 = counts[counts$variance >= quantile(counts$variance, c(.50)), ] #50% most variable genes
#counts2$variance <- NULL
#counts=counts2

x$filtr=NULL
#outlier removal if necessary

metaData <- read.csv('a.saxatilis_attributes.csv', header = TRUE, sep = ",") 
head(metaData)
nrow(metaData)


totalCounts=colSums(x) 
min(totalCounts) #910180
max(totalCounts) #2033474
mean(totalCounts) #1588113

dds <- DESeqDataSetFromMatrix(
  countData = x, # Our prepped data frame with counts
  colData = metaData, # Data frame with annotation for our samples
  design = ~1 # Here we are not specifying a model
)
dds<-dds[rowSums(counts(dds))>0,]
summary(dds) #10658


#to obtain variance stabilized data follow: 
head(dds)
vsd=vst(dds)
head(vsd)

# Check that the data has the correct format for many functions operating on multiple sets:
#WGCNA THE DATA NEEDS TO BE TRANSPOSED to make dendrogram

gsg = goodSamplesGenes(x, verbose = 3);
gsg$allOK

xx <- assay(vsd) %>%
  t() # Transpose this data

dim(vsd) #10658    40
dim(xx) #40 10658
rownames(xx) #make sure the rownames are actually 

sampleTree = hclust(dist(xx), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")


traitData = read.csv("a.saxatilis_traits.csv",header=TRUE);
dim(traitData)
names(traitData)


fishID = rownames(xx);
traitRows = match(fishID, traitData$FishID);
datTraits = traitData[traitRows, -1];
rownames(datTraits) = traitData[traitRows, 1];
collectGarbage();

#store the data in an R function 

save(xx, datTraits, file = "ASax-wgcna-step1.RData")

##To Find blocks of modules 8, use the data that is not transposed 
# Choose a set of soft-thresholding powers
lnames= load(file="ASax-wgcna-step1.RData")
lnames
exprSize = checkSets(xx, checkStructure = TRUE);
nSets = exprSize$nSets;
nSets = checkSets(xx, checkStructure = TRUE)$nSets

powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(xx, powerVector = powers, verbose = 6)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#the curve is initially stabilized at 10 - with xx dataframe
#build networks
net = blockwiseModules(xx, power = 10, corType="bicor", networkType = "signed",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = FALSE,
                       deepSplit = 2,
                       pamRespectsDendro = F,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "ASax-net",
                       verbose = 5, maxBlockSize = 5000)
names(net)




moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
consMEs = net$MEs;
consTree = net$dendrograms[[1]]

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

save(consMEs, moduleLabels, moduleColors, consTree, file = "ASax-NetworkConstruction-auto.RData")
######## 
#Clusters to All Samples (not very clear pattern)
module_df <- data.frame(gene_id = names(net$colors),
                        colors = labels2colors(net$colors))
module_df
module_df[1:5,]
write.table(module_df, file ="ASax-gene_modules.txt",sep = "\t", quote=FALSE)

write.csv(module_df, file="gene_modules.csv", quote = FALSE)

MEs0 <- moduleEigengenes(xx, mergedColors)$eigengenes
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)
MEs0$treatment = row.names(MEs0)

mME = MEs0 %>% pivot_longer(-treatment) %>% mutate(
  name = gsub("ME", "", name),
  name = factor(name, levels = module_order))

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

########
#Cluster Dendrogram Dissimilarity and Merged Dynamic

softPower = 10; #10 with outliers included
adjacency = adjacency(xx, power = softPower, type='signed');
TOM = TOMsimilarity(adjacency, TOMType='signed');
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
#pdf(file=Dendrogram_signed_BM10.pdf, width=20, height=20)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
#dev.off()
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
#dynamicMods


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
#dynamicColors

# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
#pdf(file=Dendrogram_signed_BM10_colors.pdf, width=20, height=20)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(xx, colors = dynamicColors)
MEList$eigengenes #gives you the eigenes by sample 
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
#pdf(file=ClusteringEigengenes.pdf, width=20, height=20)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#dev.off()

MEDissThres = 0.25 
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(xx, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors; #Change "net" to "merge" to do the rest with the merged colors!!!!
table(mergedColors)
mergedColors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#sizeGrWindow(12, 9)
#pdf(file = "DendroAndColors_sft6_bm10.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# MERGING: Rename to moduleColors that are similar
moduleColors = mergedColors
table(moduleColors)
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "ASax-RNA_WGCNA_networkConstruct_signed.RData")
# Define numbers of genes and samples
nGenes = ncol(xx);
nSamples = nrow(xx);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(xx, moduleColors)$eigengenes 
MEs = orderMEs(MEs0)

######################
#correlations of traits with eigengenes

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitCor 

## Will display correlations and their p-values
#textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
textMatrix =  paste(signif(moduleTraitPvalue, 1))
dim(textMatrix) = dim(moduleTraitCor)
par(mfrow = c(1,1))
par(mar=c(1,1,1,1))
# Display the correlation values within a heatmap plot
pdf(file="trait-cor-heatmap.pdf", width=7, height=7)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               xLabelsAngle = 0,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.75,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

#################################################################################
### make heatmap of specific modules of interest
head(MEs)

#Which modules have biggest differences across treatment groups?
#First we double check that our samples are still in order.
all.equal(metaData$FishID, rownames(MEs))

# Create the design matrix from the `treatment` variable
des_mat <- model.matrix(~ metaData$Treatment)

#Run linear model on each module. 
#Limma wants our tests to be per row, so we also need to transpose so the eigengenes are rows
# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(MEs), design = des_mat)
# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(MEs)) %>%
  tibble::rownames_to_column("module")

head(stats_df)


#retrieve gene list for the modules, can be modified in excel to run GO analyses
annot=read.table(file="GO_table.tab", header=T, sep='\t')
probes = colnames(xx)
probes2annot = match(probes,annot$gene)
head(probes2annot)

datGS.Traits=data.frame(cor(xx,datTraits,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datME=moduleEigengenes(xx,moduleColors)$eigengenes
datKME=signedKME(xx, datME, outputColumnName="MM.")
datOutput=list(gene=colnames(xx),colors = moduleColors)
datOutput=as.vector(datOutput)
datOutput$gene = as.character(datOutput$gene)

datOutput0 = as.data.frame(datOutput)

#As a sanity check, let's use ggplot to see what module XX's eigengene looks like between treatment groups.
#First we need to set up the module eigengene for this module with the sample metadata labels we need.
module_darkgrey_df <- MEs %>%
  tibble::rownames_to_column("FishName") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metaData %>%
                      dplyr::select(FishID, Treatment),
                    by = c("FishName" = "FishID")
  )

ggplot(
  module_darkgrey_df,
  aes(
    x = Treatment,
    y = MEdarkgrey,
    color = Treatment
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()




#If you want to know which of your genes make up a modules, you can look at the $colors slot.
gene_module_key <- as.data.frame(datOutput) %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(colors = paste0("ME", colors))
gene_module_key <- gene_module_key %>%
  dplyr::filter(colors == "MEdarkgrey")

head(gene_module_key)

#Save the gene to module key ato a TSV file
#this list was used to identify annotated genes on the command line
readr::write_csv(gene_module_key, "brain_wgcna_module_darkgrey.csv")

write.csv(net$colors, "module-genes_brain_WGCNA.csv")

#############################################################################

#Make a custom heatmap function for specific modules

make_module_heatmap <- function(module_name,
                                expression_mat = xx,
                                metadata_df = metaData,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME0"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with Treatment and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its Treatment
  #module_eigengene <- module_eigengenes_df %>%
  #dplyr::select(all_of(module_name)) %>%
  #tibble::rownames_to_column("FishID")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(Treatment, FishID) %>%
    # Add on the eigengene expression by joining with sample IDs
    #dplyr::inner_join(module_eigengene, by = "FishID") %>%
    # Arrange by patient and time point
    dplyr::arrange(Treatment, FishID) %>%
    # Store sample
    tibble::column_to_rownames("FishID")
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    treatment = col_annot_df$Treatment,
    # Add annotation barplot
    #module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in treatment
    col = list(treatment = c("CHC" = "#009E73", "CLC" = "#56B4E9", "HWHC" = "#E69F00", "HWLC" = "#D55E00"))
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(colors == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("blue", "white", "red")
  )
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = FALSE,
                                     show_column_names = FALSE
  )
  
  # Return heatmap
  return(heatmap)
}

#Specify which module you want to investigate
module_darkgrey_heatmap <- make_module_heatmap(module_name = "MEdarkgrey")
# Print out the plot
module_darkgrey_heatmap

#Save plot to png
pdf("darkgrey_heatmap.pdf")
module_darkgrey_heatmap
dev.off()


