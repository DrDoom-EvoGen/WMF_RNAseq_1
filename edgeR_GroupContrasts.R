# Moving quantification files into edgeR and setting up edgeR files and contrasts

## This version uses each grouped contols in contrasts

library(edgeR)

# Import quantification file from TxImport step 
txi <- readRDS("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/TxImport.RDS")

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
group <- c('E_C','E_C','E_C','E_C','E_C','E_C',
           'E_C','E_C','E_C','E_C','E_C','E_C',
           'E_T_12','E_T_12','E_T_12','E_T_24','E_T_24','E_T_24',
           'E_T_48','E_T_48','E_T_48','E_T_96','E_T_96','E_T_96',
           'H_C','H_C','H_C','H_C','H_C','H_C',
           'H_C','H_C','H_C','H_C','H_C','H_C',
           'H_T_12','H_T_12','H_T_12','H_T_24','H_T_24','H_T_24',
           'H_T_48','H_T_48','H_T_48','H_T_96','H_T_96','H_T_96')

y <- DGEList(counts = cts, group = group)
y <- scaleOffset(y, normMat)

head(y)

# filtering
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

plotMDS(y)

# for plotting other axes of MDS
#MDS <- plotMDS(y,ndim=5,plot=FALSE)
#dims <- MDS$cmdscale.out

# y is now ready for estimate dispersion functions see edgeR User's Guide

design <- model.matrix(~0+group, data = y$samples)

colnames(design) <- unique(group)

y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)

# Now we build the contrasts we want to make. These can be informed by the MDS plot

my.contrasts <- makeContrasts(
  control = E_C-H_C
  E_12 = E_T_12-E_C,
  E_24 = E_T_24-E_C,
  E_48 = E_T_48-E_C,
  E_96 = E_T_96-E_C,
  H_12 = H_T_12-H_C,
  H_24 = H_T_24-H_C,
  H_48 = H_T_48-H_C,
  H_96 = H_T_96-H_C,
  GENO_12 = (E_T_12-E_C)-(H_T_12-H_C),
  GENO_24 = (E_T_24-E_C)-(H_T_24-H_C),
  GENO_48 = (E_T_48-E_C)-(H_T_48-H_C),
  GENO_96 = (E_T_96-E_C)-(H_T_96-H_C),
  levels=design)

saveRDS(fit, file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/Group_fit.RDS")
saveRDS(my.contrasts, file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/GroupContrasts.RDS")


