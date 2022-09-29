# WMF_RNAseq_1
Road map and documentation of RNA-seq analysis of one hybrid (from Hayden Lake, Idaho) and one Eurasian Watermilfoil (from Couer D'Alene Lake, Idaho) genotype 

Initial steps were run on a Linux-based high performance computing cluster at Montana State University. This cluster uses a SLURM job management system. Files run this way have ".sbatch" file extensions. The rest of the files are R scripts.

"bbadaptertrim.sbatch" = Remove adapters from raw reads with BBduk.

"bbqualitytrim.sbatch" = Trim poor quality reads with BBduk.

"fastqc.sbatch" = Quality check of cleaned reads.

"salmonindex.sbatch" = Index the Trinity assembled transcriptome (see WMF_DeNovoTranscriptome) for Salmon to align experimental too.

"salmonalign.sbatch" = Align experimental reads to Trinity reference with Salmon. The output from this alignment are quantification files for each experimental sample used in this study.

The following analyses were done in R.

"TxImport.R" = An R script that moves Salmon quantification files into R and aggregates them into gene counts using TXimport package.

"edgeR_Contrasts.R" = An R script that uses edgeR to calculate a full model and set contrasts for differential expression of genes. 

"GoseEnrichment.R" = An R script for enrichment analysis using the GOseq package and to determine shared differentially expressed genes and plot the numbers in venn diagrams with the VennDetail package.

DEcompare.R =  An R script that aggregates all log-fold-change files and constructs heatmaps of Gene Ontology terms of interest using the PHeatmaps package.

DoseResponse.R = An R script for calculating dose response models using the drc package.

Manuscript published at: 