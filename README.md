#### This repository corresponds with transcriptome annotation and data analysis for the manuscript "Molecular plasticity to ocean warming and habitat loss in a coral reef fish"
- linked doi when available

![Treatments_ASax](https://github.com/allyswank/Abudefduf-Brains/assets/91483379/30651a8f-9bd3-4c5b-ac8b-d9c7d1580ccc)


# Contents
* [Tag-Seq Processing](#Sequence_Processing)
* [Annotation](#Annotation)
* [Differential Gene Expression](#Differential_Gene_Expression)
* [WGCNA](#WGCNA)
* [ELISA](#ELISA)

---

## Sequence Processing
Total RNA samples from *Abudefduf saxatilis* were shipped to the University of Texas at Austin Genomic Sequencing and Analysis Facility for Tag-Seq library preparation (following Meyer et al. 2011) and sequencing using the Illumina NovaSeq6000. 

* Raw Tag-seq reads and transcriptome assembly accessions can be found under SRA Bioproject PRJNA1009081.
* Sequence quality was assessed with FastQC v0.11.9: `fastqc` and `fastqc_trimmed`
* We used the program cutadapt v1.13 to trim Illumina adapters: `clean`
* Reads were mapped to the reference transcriptome with Bowtie2 v2.3.: `maps`

Custom scripts and a comprehensive review of the read-sequencing pipeline are available at: https://github.com/z0on/tag-based_RNAseq



## Annotation
The assembled transcriptome was annotated using the program BLAST+ ver 2.11.0 (Camacho et al. 2009) to search for sequence similarity against closely related gene sets on Ensembl (Cunningham et al. 2022). After a BLAST search, matches were filtered for best matches to each transcript, resulting in a file of remaining unnanotated transcripts and a complete annotation table (`annotated_transcripts.txt`).

We repeated BLAST searches and filtering steps iteratively accross Ensembl gene sets in the following order for transcripts that remained unannotated: 
* *Amphiprion percula* (orange clownfish): `BLAST_Nemo.sh` & `filter_Nemo.sh`
* *Acanthochromis polyacanthus* (spiny chromis): `BLAST_Chromis.sh` & `filter_Chromis.sh`
* *Danio rerio* (zebrafish): `BLAST_Zebra.sh` & `filter_Zebra.sh`
* *Stegastes partitus* (bicolor damselfish): `BLAST_Steg.sh` & `filter_Steg.sh`
* *Amphiprion ocellaris* (clown anemonefish): `BLAST_Clown.sh` & `filter_Clown.sh`

Any remaining unannotated contigs were compared to the UniProt Swiss-Prot protein database (Apweiler et al. 2004): `BLAST_uniprot.sh` & `filter_Uniprot.sh`


## Differential Gene Expression
The analyses of differential gene expression were conducted with DESeq2 ver 1.34.0 (Love et al. 2014) in RStudio ver 4.1.3 (RStudio Team, 2019) using the `Asaxatilis_RNAseq.R` script. 

Estimates GO category enrichment were identified using Log2Fold values of all gene counts using `GO_MWU.R`. 

All scripts were adapted from https://github.com/z0on/GO_MWU 

## WGCNA
A weighted gene co-expression network analysis (WGCNA) was used to identify clusters of highly correlated genes within the experimental conditions (Langfelder and Horvath 2008) using the `WGCNA.R` script. 

GO enrichment estimates were made on all significant gene modules by adapting `GO_MWU.R`. 


## ELISA
A protein carbonyl enzyme-linked immunosorbent assay (ELISA) kit was used to quantify the protein oxidation levels of liver and muscle tissues separately. Because they are well conserved and easily detected, protein carbonyl formations provide insight to oxidative stress levels across tissues and taxa (Birnie-Gauvin et al. 2017).

All statistical analyses for raw protein carbonyl data (`A.saxatilis_ELISA-w-attributes.csv`) are in the script `Asaxatilis_ELISA.R`.
