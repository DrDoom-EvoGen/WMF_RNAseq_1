# Moving quantification files into edgeR and setting up edgeR files and contrasts

## This version uses each grouped contols in contrasts

library(edgeR)
library(ggforce)

# Import quantification file from TxImport step 
txi <- readRDS("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/TxImport.RDS")

names(txi)
head(txi$counts)
head(txi$length)

# Use Counts and Lengths to build DGE object
cts <- txi$counts
normMat <- txi$length

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.

eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
group <- c('E_C_12','E_C_12','E_C_12','E_C_24','E_C_24','E_C_24',
           'E_C_48','E_C_48','E_C_48','E_C_96','E_C_96','E_C_96',
           'E_T_12','E_T_12','E_T_12','E_T_24','E_T_24','E_T_24',
           'E_T_48','E_T_48','E_T_48','E_T_96','E_T_96','E_T_96',
           'H_C_12','H_C_12','H_C_12','H_C_24','H_C_24','H_C_24',
           'H_C_48','H_C_48','H_C_48','H_C_96','H_C_96','H_C_96',
           'H_T_12','H_T_12','H_T_12','H_T_24','H_T_24','H_T_24',
           'H_T_48','H_T_48','H_T_48','H_T_96','H_T_96','H_T_96')

y <- DGEList(counts = cts, group = group)
y <- scaleOffset(y, normMat)

head(y)

# filtering
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/MDSplot.pdf",
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
plotMDS(y, col = as.numeric(as.factor(group)))
dev.off()

# for plotting other axes of MDS
MDS <- plotMDS(y,ndim=5,plot=FALSE, top = 2000)
dims <- data.frame(MDS$cmdscale.out)
dims$label <- row.names(dims)
dims[c('Genotype', 'Treatment', 'Time', 'Rep')] <- str_split_fixed(dims$label , '_', 4)

E_dims <- dims %>%
  filter(Genotype == "E")

H_dims <- dims %>%
  filter(Genotype == "H")

ggplot(dims, aes(x = X1, y = X2)) +
  geom_point(data = E_dims, aes(shape = Treatment), color = "grey", cex = 4) +
  geom_point(data = H_dims, aes(shape = Treatment), color = "black", cex = 4) +
  theme_minimal() + 
  xlab("Leading log-FC Dimension 1") + 
  ylab("Leading log-FC Dimension 2") +
  geom_mark_ellipse(data = E_dims, aes(x = X1, y = X2, color = Time, linetype = Treatment)) +
  geom_mark_ellipse(data = H_dims, aes(x = X1, y = X2, color = Time, linetype = Treatment))

ggsave("MDS.Plot.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  width = 5, 
  height = 5,
  dpi = 600)

# y is now ready for estimate dispersion functions see edgeR User's Guide

design <- model.matrix(~0+group, data = y$samples)

colnames(design) <- unique(group)

y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)

# Now we build the contrasts we want to make. These can be informed by the MDS plot

my.contrasts <- makeContrasts(
  control = (E_C_12 + E_C_24 + E_C_48 + E_C_96)/4 - (H_C_12 + H_C_24 + H_C_48 + H_C_96)/4,
  E_12 = E_T_12 - E_C_12,
  E_24 = E_T_24 - E_C_24,
  E_48 = E_T_48 - E_C_48,
  E_96 = E_T_96 - E_C_96,
  H_12 = H_T_12 - H_C_12,
  H_24 = H_T_24 - H_C_24,
  H_48 = H_T_48 - H_C_48,
  H_96 = H_T_96 - H_C_96,
  HAT_12 = (E_T_12 - E_C_12) - (H_T_12 - H_C_12),
  HAT_24 = (E_T_24 - E_C_24) - (H_T_24 - H_C_24),
  HAT_48 = (E_T_48 - E_C_48) - (H_T_48 - H_C_48),
  HAT_96 = (E_T_96 - E_C_96) - (H_T_96 - H_C_96),
  E_treat = (E_T_12 + E_T_24 + E_T_48 + E_T_96)/4 - (E_C_12 + E_C_24 + E_C_48 + E_C_96)/4,
  H_treat = (H_T_12 + H_T_24 + H_T_48 + H_T_96)/4 - (H_C_12 + H_C_24 + H_C_48 + H_C_96)/4,
  levels=design)

saveRDS(fit, file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/Group_fit.RDS")
saveRDS(my.contrasts, file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/GroupContrasts.RDS")
