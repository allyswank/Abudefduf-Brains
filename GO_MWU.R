
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 

# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu

################################################################
# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.
setwd("C:/Users/allys/Box/Bernal_lab/Ally/Chapter1/Analyses/GO")


#install.packages("ape")
library("ape")



# Edit these to match your data file names: 
input="WGCNA_grey60_fishers.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="GO_table.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="CC" # either MF, or BP, or CC
source("gomwu.functions.R")



# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath="C:/Strawberry/perl/bin/perl.exe",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25,
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.


# ----------- Plotting results

results=gomwuPlot(input,goAnnotations,goDivision,
                  #absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)

results

 # manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]


# ------- extracting representative GOs

# this module chooses GO terms that best represent *independent* groups of significant GO terms

pcut=1e-2 # adjusted pvalue cutoff for representative GO
hcut=0.9 # height at which cut the GO terms tree to get "independent groups". 

# plotting the GO tree with the cut level (un-remark the next two lines to plot)
plot(results[[2]],cex=0.6)
abline(h=hcut,col="red")

# cutting
ct=cutree(results[[2]],h=hcut)
annots=c();ci=1
for (ci in unique(ct)) {
  message(ci)
  rn=names(ct)[ct==ci]
  obs=grep("obsolete",rn)
  if(length(obs)>0) { rn=rn[-obs] }
  if (length(rn)==0) {next}
  rr=results[[1]][rn,]
  bestrr=rr[which(rr$pval==min(rr$pval)),]
  best=1
  if(nrow(bestrr)>1) {
    nns=sub(" .+","",row.names(bestrr))
    fr=c()
    for (i in 1:length(nns)) { fr=c(fr,eval(parse(text=nns[i]))) }
    best=which(fr==max(fr))
  }
  if (bestrr$pval[best]<=pcut) { annots=c(annots,sub("\\d+\\/\\d+ ","",row.names(bestrr)[best]))}
}

mwus=read.table(paste("MWU",goDivision,input,sep="_"),header=T)
bestGOs=mwus[mwus$name %in% annots,]
bestGOs

write.csv(mwus, file="dr_BP_LC_CvsHW.csv", quote=FALSE)

################################################################################
#Formatting files from WGCNA for Fisher's Test

#genes in LC_CvsHW

all_genes = read.csv("counted.txt") #all genes in count table from raw RNA-seq data
Wald_LC_CvsHW = read.csv("Wald_LC_CvsHW.csv") #genes in pairwise comparison that were exclusive to low complexity (HC-CvsHW : LC-CvsHW)

#create binary score column for fisher's test
all_genes$score <- all_genes$gene %in% Wald_LC_CvsHW$gene
all_genes$score[all_genes$score=="TRUE"]<-"1"
all_genes$score[all_genes$score=="FALSE"]<-"0"

write.csv(all_genes, file="Fishers_LC_CvsHW.csv", quote=FALSE)

#WGCNA modules

all_genes_lightcyan = read.csv("counted.txt") #all genes in count table from raw RNA-seq data
lightcyan_genes = read.csv("brain_wgcna_module_lightcyan.tsv", sep = "\t") #genes in lightcyan WGCNA module
lightcyan_genes <- as.data.frame(lightcyan_genes[,1])
colnames(lightcyan_genes) <- c("gene")


#create binary score column for fisher's test
all_genes_lightcyan$score <- all_genes_lightcyan$gene %in% lightcyan_genes$gene
all_genes_lightcyan$score[all_genes_lightcyan$score=="TRUE"]<-"1"
all_genes_lightcyan$score[all_genes_lightcyan$score=="FALSE"]<-"0"

length(which(all_genes_lightcyan$score=="1")) #dummy check
write.csv(all_genes_lightcyan, file="lightcyan_WGCNA_fishers.csv", quote=FALSE)



