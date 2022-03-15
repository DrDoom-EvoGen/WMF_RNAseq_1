# WMF_RNAseq_1
RNA-seq analysis of 1 hybrid (Hayden) and 1 EWM (CDA E)

1. Remove adapters from raw reads with BBduk
2. Trim poor quality reads with BBduk
3. Align experimental reads to Trinity reference with Salmon
4. Import quantification files from Salmon into R using TXimport
5. Make DGE in edgeR 
6. Analyze DE using edgeR
7. Enrichment analysis using GOseq
8. Venn diagrams with VennDetail
9. Heatmaps


Published at:  