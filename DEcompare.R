
## Compare genes responding to the Treatment at each timepoint and between genotypes:

# Load in R packages
library(edgeR)
library(readr)
library(topGO)
library(dplyr)
library(goseq)
library(ggplot2)
library(tidyverse)
library(data.table)
library(readr)
library(stringr)
library(pheatmap)
library(randomcoloR)
library(VennDetail)

# Load in data files from previous scripts:
Control_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_control.csv")

E12_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_E12.csv")
E24_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_E24.csv")
E48_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_E48.csv")
E96_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_E96.csv")

H12_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_H12.csv")
H24_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_H24.csv")
H48_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_H48.csv")
H96_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_H96.csv")

HAT12_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_HAT12.csv")
HAT24_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_HAT24.csv")
HAT48_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_HAT48.csv")
HAT96_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_HAT96.csv")

E_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_E_treat.csv")
H_DE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_H_treat.csv")

txi <- readRDS("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/TxImport.RDS")

shared12 <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_12.csv")
shared24 <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_24.csv")
shared48 <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_48.csv")
shared96 <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_96.csv")
sharedE <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_E.csv")
sharedH <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_H.csv")
sharedEarly <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_Early.csv")
sharedTreat <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_Treat.csv")

Control.sig <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/Control_sig.csv")

E12.sig <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/E12_sig.csv")
E24.sig <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/E24_sig.csv")
E48.sig <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/E48_sig.csv")
E96.sig <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/E96_sig.csv")

H12.sig <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/H12_sig.csv")
H24.sig <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/H24_sig.csv")
H48.sig <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/H48_sig.csv")
H96.sig <- read.csv("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/H96_sig.csv")


# Filter for significantly (BH corrected P < or = 0.05) DEGs
Control_DE.sig <- subset(Control_DE, X %in% Control.sig$X)

E12_DE.sig <- subset(E12_DE, X %in% E12.sig$X)
E24_DE.sig <- subset(E24_DE, X %in% E24.sig$X)
E48_DE.sig <- subset(E48_DE, X %in% E48.sig$X)
E96_DE.sig <- subset(E96_DE, X %in% E96.sig$X)

H12_DE.sig <- subset(H12_DE, X %in% H12.sig$X)
H24_DE.sig <- subset(H24_DE, X %in% H24.sig$X)
H48_DE.sig <- subset(H48_DE, X %in% H48.sig$X)
H96_DE.sig <- subset(H96_DE, X %in% H96.sig$X)


# Combine all significant DE on the Trinity Gene
# ** E12 is not appended in this table
all_DE <- E12_DE.sig %>%
full_join(E24_DE.sig, by = 'X', copy = TRUE, suffix = c("",".E24")) %>%
full_join(E48_DE.sig, by = 'X', copy = TRUE, suffix = c("",".E48")) %>%
full_join(E96_DE.sig, by = 'X', copy = TRUE, suffix = c("",".E96")) %>%
full_join(H12_DE.sig, by = 'X', copy = TRUE, suffix = c("",".H12")) %>%
full_join(H24_DE.sig, by = 'X', copy = TRUE, suffix = c("",".H24")) %>%
full_join(H48_DE.sig, by = 'X', copy = TRUE, suffix = c("",".H48")) %>%
full_join(H96_DE.sig, by = 'X', copy = TRUE, suffix = c("",".H96")) %>%
full_join(Control_DE.sig, by = 'X', copy = TRUE, suffix = c("",".Control"))

all_DE <- all_DE %>%
unique()

ano <- read.csv("/Users/Greg/Documents/GitHub/WMF_DeNovoTranscriptome/results/TrinotateTranscriptome/Trinotate.csv")

blastx.ano <- subset(ano, select = c(X.gene_id, sprot_Top_BLASTX_hit))
blastx.sprot <- as.data.frame(str_split_fixed(blastx.ano$sprot_Top_BLASTX_hit, "_", 2))
blastx <- data.frame(matrix(ncol = 2, nrow = 579697))
colnames(blastx) <- c('gene', 'sprot')
blastx$gene <- blastx.ano$X.gene_id
blastx$sprot <- blastx.sprot$V1
blastx <- blastx %>% 
unique()
blastx$sprot[blastx$sprot=="."] <- NA
blastx <- na.omit(blastx)



## Heat maps using GO enriched terms

go_data <- topGO::readMappings("/Users/Greg/Documents/GitHub/WMF_DeNovoTranscriptome/results/TrinotateTranscriptome/gene_ontology_proc.txt")
n_go <- lapply(go_data, length) %>% unlist
gene_go <- data.frame(gene_id=rep(names(go_data), n_go),
                      go_term=unlist(go_data)) %>%
  arrange(gene_id) %>%
  unique()
head(gene_go)

logFC <- subset(all_DE, select = c(X, logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control))
write.csv(logFC, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/logFC.csv")

my.breaks <- c(seq(-12, -.01, by=0.1), seq(0, 12, by=0.1)) 
my.colors <- c(colorRampPalette(colors = c("dark blue", "blue", "light blue", "white"))(length(my.breaks)/2), colorRampPalette(colors = c("white", "yellow", "orange", "red"))(length(my.breaks)/2))

annotation_col = data.frame(
                    Genotype = factor(c("EWM", "EWM", "EWM", "EWM", "hybrid", "hybrid", "hybrid", "hybrid", "EWM - hybrid")))

rownames(annotation_col) = c("logFC", "logFC.E24", "logFC.E48", "logFC.E96", "logFC.H12", "logFC.H24", "logFC.H48", "logFC.H96", "logFC.Control")

ann_colors = list(
    Genotype = c(EWM = "white", hybrid = "grey", 'EWM - hybrid' = "black"))

### Main Heatmaps

## 9-cis-epoxycarotenoid dioxygenase activity
nced <- gene_go %>% filter(gene_go$go_term == 'GO:0045549')
nced_logFC <- subset(logFC, X %in% nced$gene_id)
nced.mtx <- as.matrix(subset(nced_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(nced.mtx) <- nced_logFC$X

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/heats/NCED.pdf",   # Phaetmaps wont size this plot right with the "filename" option
    width = 4, # The width of the plot in inches
    height = 6) # The height of the plot in inches

pheatmap(nced.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 6, 
         fontsize_row = 5, 
         fontsize = 8,
         show_colnames = TRUE,
         show_rownames = FALSE,
         legend = TRUE,
         legend_labels = c("Log-Fold_Change"),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         annotation_names_col = FALSE,
         labels_col = c("12", "24", "48", "96", "12", "24", "48", "96", "ALL"),
         main = "NCED Activity",
         border_color=NA)
dev.off()

## Abscisic Acid Signaling Pathway
aba2 <- gene_go %>% filter(gene_go$go_term == 'GO:0009738')
aba2_logFC <- subset(logFC, X %in% aba2$gene_id)
aba2.mtx <- as.matrix(subset(aba2_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(aba2.mtx) <- aba2_logFC$X

pheatmap(aba2.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 2, 
         fontsize_row = 5, 
         fontsize = 8,
         show_colnames = TRUE,
         show_rownames = FALSE,
         legend = TRUE,
         legend_labels = c("Log-Fold_Change"),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         annotation_names_col = FALSE,
         labels_col = c("12", "24", "48", "96", "12", "24", "48", "96", "ALL"),
         main = "Abscisic Acid Signaling Pathway",
         filename = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/heats/ABASignaling.pdf")


## Photosynthesis
photo <- gene_go %>% filter(gene_go$go_term == 'GO:0015979')
photo_logFC <- subset(logFC, X %in% photo$gene_id)
photo.mtx <- as.matrix(subset(photo_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(photo.mtx) <- photo_logFC$X

pheatmap(photo.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 2, 
         fontsize_row = 5, 
         fontsize = 8,
         show_colnames = TRUE,
         show_rownames = FALSE,
         legend = TRUE,
         legend_labels = c("Log-Fold_Change"),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         annotation_names_col = FALSE,
         labels_col = c("12", "24", "48", "96", "12", "24", "48", "96", "ALL"),
         main = "Photosynthesis",
         filename = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/heats/Photosynthesis.pdf")


## Ethylene Activated Signaling 
ethylene <- gene_go %>% filter(gene_go$go_term == 'GO:0009873')
ethylene_logFC <- subset(logFC, X %in% ethylene$gene_id)
ethylene.mtx <- as.matrix(subset(ethylene_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(ethylene.mtx) <- ethylene_logFC$X

pheatmap(ethylene.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 2, 
         fontsize_row = 5, 
         fontsize = 8,
         show_colnames = TRUE,
         show_rownames = FALSE,
         legend = TRUE,
         legend_labels = c("Log-Fold_Change"),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         annotation_names_col = FALSE,
         labels_col = c("12", "24", "48", "96", "12", "24", "48", "96", "ALL"),
         main = "Ethylene-Activated Signaling",
         filename = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/heats/Ethylene.pdf")









### Additional Heatmaps not in manuscript

## Auxin efflux
AuxEfflux <- gene_go %>% filter(gene_go$go_term == 'GO:0010315')
AuxEfflux_logFC <- subset(logFC, X %in% AuxEfflux$gene_id)
AuxEfflux.mtx <- as.matrix(subset(AuxEfflux_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(AuxEfflux.mtx) <- AuxEfflux_logFC$X

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/AuxinEfflux.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 12) # The height of the plot in inches

pheatmap(AuxEfflux.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 3, 
         fontsize_row = 3, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "Auxin Efflux GO:0010315")
dev.off()


## Auxin Transport
AuxTrans <- gene_go %>% filter(gene_go$go_term == 'GO:0060918')
AuxTrans_logFC <- subset(logFC, X %in% AuxTrans$gene_id)
AuxTrans.mtx <- as.matrix(subset(AuxTrans_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(AuxTrans.mtx) <- AuxTrans_logFC$X

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/AuxinTransport.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 12) # The height of the plot in inches

pheatmap(AuxTrans.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 3, 
         fontsize_row = 3, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "Auxin Transport GO:0060918")
dev.off()

## Hormone Transport
HormoneTrans <- gene_go %>% filter(gene_go$go_term == 'GO:0009914')
HormoneTrans_logFC <- subset(logFC, X %in% HormoneTrans$gene_id)
HormoneTrans.mtx <- as.matrix(subset(HormoneTrans_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(HormoneTrans.mtx) <- HormoneTrans_logFC$X

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/HormoneTransport.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 12) # The height of the plot in inches

pheatmap(HormoneTrans.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 3, 
         fontsize_row = 3, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "Hormone Transport GO:0009914")
dev.off()

## regulation of auxin mediated signaling pathway
AuxSigReg <- gene_go %>% filter(gene_go$go_term == 'GO:0010928')
AuxSigReg_logFC <- subset(logFC, X %in% AuxSigReg$gene_id)
AuxSigReg.mtx <- as.matrix(subset(AuxSigReg_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(AuxSigReg.mtx) <- AuxSigReg_logFC$X

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/AuxinSigRegulation.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 12) # The height of the plot in inches

pheatmap(AuxSigReg.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 3, 
         fontsize_row = 3, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "regulation of auxin mediated signaling pathway GO:0010928")
dev.off()

## auxin transmembrane transporter activity
AuxTransTrans <- gene_go %>% filter(gene_go$go_term == 'GO:0080161')
AuxTransTrans_logFC <- subset(logFC, X %in% AuxTransTrans$gene_id)
AuxTransTrans.mtx <- as.matrix(subset(AuxTransTrans_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(AuxTransTrans.mtx) <- AuxTransTrans_logFC$X

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/AuxTransTransporter.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 12) # The height of the plot in inches

pheatmap(AuxTransTrans.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 3, 
         fontsize_row = 3, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "auxin transmembrane transporter activity GO:0080161")
dev.off()

## auxin efflux transmembrane transporter activity
AuxEffluxTrans <- gene_go %>% filter(gene_go$go_term == 'GO:0010329')
AuxEffluxTrans_logFC <- subset(logFC, X %in% AuxEffluxTrans$gene_id)
AuxEffluxTrans.mtx <- as.matrix(subset(AuxEffluxTrans_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(AuxEffluxTrans.mtx) <- AuxEffluxTrans_logFC$X

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/AuxEffluxTransporter.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 12) # The height of the plot in inches

pheatmap(AuxEffluxTrans.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 3, 
         fontsize_row = 3, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "auxin efflux transmembrane transporter activity GO:0010329")
dev.off()

##Auxin activated signalling
HormoneSig <- gene_go %>% filter(gene_go$go_term == 'GO:0009734')
HormoneSig_logFC <- subset(logFC, X %in% HormoneSig$gene_id)
HormoneSig.mtx <- as.matrix(subset(HormoneSig_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(HormoneSig.mtx) <- HormoneSig_logFC$X

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/AuxinSignaling.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 12) # The height of the plot in inches

pheatmap(HormoneSig.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 1.5, 
         fontsize_row = 3, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "Auxin Activated Signalling")
dev.off()



##Auxin homeostasis
AuxStasis <- gene_go %>% filter(gene_go$go_term == 'GO:0010252')
AuxStasis_logFC <- subset(logFC, X %in% AuxStasis$gene_id)
AuxStasis.mtx <- as.matrix(subset(AuxStasis_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(AuxStasis.mtx) <- AuxStasis_logFC$X

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/AuxinHomeostasis.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 12) # The height of the plot in inches

pheatmap(AuxStasis.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 5, 
         fontsize_row = 5, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "Auxin Homeostasis")
dev.off()

##ABA metabolism
ABAmetab <- gene_go %>% filter(gene_go$go_term == 'GO:0009687')
ABAmetab_logFC <- subset(logFC, X %in% ABAmetab$gene_id)
ABAmetab.mtx <- as.matrix(subset(ABAmetab_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(ABAmetab.mtx) <- ABAmetab_logFC$X

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/ABAmetab.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 12) # The height of the plot in inches

pheatmap(ABAmetab.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 5, 
         fontsize_row = 5, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col)
dev.off()

##negative regulation of abscisic acid-activated signaling pathway
ABAnegReg <- gene_go %>% filter(gene_go$go_term == 'GO:0009788')
ABAnegReg_logFC <- subset(logFC, X %in% ABAnegReg$gene_id)
ABAnegReg.mtx <- as.matrix(subset(ABAnegReg_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(ABAnegReg.mtx) <- ABAnegReg_logFC$X

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/NegativeABA.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 12) # The height of the plot in inches

pheatmap(ABAnegReg.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 5, 
         fontsize_row = 5, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "Negative Regulation of ABA Pathway")
dev.off()

## auxin transmembrane transporter activity
auxtrans <- gene_go %>% filter(gene_go$go_term == 'GO:0080161')
auxtrans_logFC <- subset(logFC, X %in% auxtrans$gene_id)
auxtrans.mtx <- as.matrix(subset(auxtrans_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(auxtrans.mtx) <- auxtrans_logFC$X
pheatmap(auxtrans.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 5, 
         fontsize_row = 5, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col)



## steroid hydroxylase activity
roid <- gene_go %>% filter(gene_go$go_term == 'GO:0008395')
roid_logFC <- subset(logFC, X %in% roid$gene_id)
roid.mtx <- as.matrix(subset(roid_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(roid.mtx) <- roid_logFC$X
pheatmap(roid.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 5, 
         fontsize_row = 5, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col)

## hormone catabolic process
hormonecatab <- gene_go %>% filter(gene_go$go_term == 'GO:0042447')
hormonecatab_logFC <- subset(logFC, X %in% hormonecatab$gene_id)
hormonecatab.mtx <- as.matrix(subset(hormonecatab_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(hormonecatab.mtx) <- hormonecatab_logFC$X
pheatmap(hormonecatab.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 5, 
         fontsize_row = 5, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col)



## Response to Abscisic Acid 
aba <- gene_go %>% filter(gene_go$go_term == 'GO:0009737')
aba_logFC <- subset(logFC, X %in% aba$gene_id)
aba.mtx <- as.matrix(subset(aba_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(aba.mtx) <- aba_logFC$X

pheatmap(aba.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 1, 
         fontsize_row = 5, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "Response to Abscisic Acid",
         filename = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/ResponseABA.pdf")

## Photosynthesis - Dark Reaction
photodark <- gene_go %>% filter(gene_go$go_term == 'GO:0019685')
photodark_logFC <- subset(logFC, X %in% photodark$gene_id)
photodark.mtx <- as.matrix(subset(photodark_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(photodark.mtx) <- photodark_logFC$X

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/PhotoDark.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 12) # The height of the plot in inches

pheatmap(photodark.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 3, 
         fontsize_row = 5, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "Photosynthesis - Dark Reaction")
dev.off()


## Response to Stress
stress <- gene_go %>% filter(gene_go$go_term == 'GO:0006950')
stress_logFC <- subset(logFC, X %in% stress$gene_id)
stress.mtx <- as.matrix(subset(stress_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(stress.mtx) <- stress_logFC$X
pheatmap(stress.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 0.5, 
         fontsize_row = 5, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "Response to Stress")


## cytokinin metabolic process# 
cytokin <- gene_go %>% filter(gene_go$go_term == 'GO:0009690')
cytokin_logFC <- subset(logFC, X %in% cytokin$gene_id)
cytokin.mtx <- as.matrix(subset(cytokin_logFC, select = c(logFC, logFC.E24, logFC.E48, logFC.E96, logFC.H12, logFC.H24, logFC.H48, logFC.H96, logFC.Control)))
rownames(cytokin.mtx) <- cytokin_logFC$X
pheatmap(cytokin.mtx, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "none", 
         breaks = my.breaks, 
         color = my.colors, 
         cellwidth = 15, 
         cellheight = 4, 
         fontsize_row = 5, 
         fontsize = 8,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "cytokinin metabolic process")


## Comparing mean control abundance between the hybrid and EWM genotype

abundance <- as.data.frame(txi$abundance)
abundance$E_C_row_mean <- rowMeans(abundance[,c(1:12)], na.rm=TRUE)
abundance$H_C_row_mean <- rowMeans(abundance[,c(25:36)], na.rm=TRUE)
abundance <- rownames_to_column(abundance, "row_name")
abundance <- subset(abundance, select = c("row_name", "E_C_row_mean", "H_C_row_mean"))
abundance$mean_sum <- abundance$E_C_row_mean + abundance$H_C_row_mean
abundance <- abundance %>% filter(abundance$mean_sum > 0)

abundance.sig <- subset(abundance, abundance$row_name %in% Control.sig$X)
abundance.sig <- subset(abundance.sig, select = c("row_name", "E_C_row_mean", "H_C_row_mean"))
write.csv(abundance.sig, "/Users/Greg/Desktop/Control_abundance_sig.csv")


ggplot(abundance, aes(E_C_row_mean, H_C_row_mean)) + 
    geom_point() +
    geom_point(data = abundance.sig, color = 'red', pch = 1) + 
    geom_abline(slope = 1, intercept = 0) +
    xlim(0,1000) +
    ylim(0,1000) + 
    xlab("EWM Mean Control Abundance") + 
    ylab("hybrid Mean Control Abundance") +
    theme_classic()

ggsave("Abundance_Corr_1000.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 600)

ggplot(abundance, aes(log(E_C_row_mean), log(H_C_row_mean))) + 
    geom_point() +
    geom_point(data = abundance.sig, color = 'red', pch = 1) + 
    geom_abline(slope = 1, intercept = 0) +
    xlim(0,12) +
    ylim(0,12) + 
    xlab("log of EWM Mean Control Abundance") + 
    ylab("log of hybrid Mean Control Abundance") +
    theme_classic()

ggsave("Abundance_Corr_log.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 600)

ggplot(abundance, aes(E_C_row_mean, H_C_row_mean)) + 
    geom_point() +
    geom_point(data = abundance.sig, color = 'red', pch = 1) + 
    geom_abline(slope = 1, intercept = 0) +
    xlim(0,250) +
    ylim(0,250) + 
    xlab("EWM Mean Control Abundance") + 
    ylab("hybrid Mean Control Abundance") +
    theme_classic()

ggsave("Abundance_Corr_250.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 600)

ggplot(abundance, aes(E_C_row_mean, H_C_row_mean)) + 
    geom_point() +
    geom_point(data = abundance.sig, color = 'red', pch = 1) + 
    geom_abline(slope = 1, intercept = 0) +
    xlim(0,50) +
    ylim(0,50) + 
    xlab("EWM Mean Control Abundance") + 
    ylab("hybrid Mean Control Abundance") +
    theme_classic()

ggsave("Abundance_Corr_50.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 600)




