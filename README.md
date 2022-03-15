# WMF_RNAseq_1
Road map and documentation of RNA-seq analysis of 1 hybrid (Hayden) and 1 Eurasian Watermilfoil  (CDA E) genotype 

Initial steps were run on a Linux-based high performance computing cluster at Montana State University. This cluster uses a SLURM job management system. Files run this way have ".sbatch" file extensions.

"bbadaptertrim.sbatch" = Remove adapters from raw reads with BBduk.

"bbqualitytrim.sbatch" = Trim poor quality reads with BBduk.

"fastqc.sbatch" = Quality check of cleaned reads.

"salmonindex.sbatch" = Index a Trinity assembled transcriptome for Salmon to align experimental too.

"salmonalign.sbatch" = Align experimental reads to Trinity reference with Salmon. The output from this alignment are quantification files for each experimental sample used in this study.

The following analyses were don in R.

4. Import quantification files from Salmon into R using TXimport
5. Make DGE in edgeR 
6. Analyze DE using edgeR
7. Enrichment analysis using GOseq
8. Venn diagrams with VennDetail
9. Heatmaps


Published at:  