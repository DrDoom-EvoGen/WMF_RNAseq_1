# Importing data from Salmon quantification files into R

library(readr)
library(tximportData)
library(tximport)

dir <- ("/Users/Greg/Desktop/quants")
list.files(dir)

samples <- (list.files(dir))
samples

files <- file.path(dir, samples, "quant.sf")
all(file.exists(files))
names(files) <- samples

tx2gene <- read_csv("/Users/Greg/Desktop/tx2gene.csv")
head(tx2gene)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

saveRDS(txi, file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/TxImport.RDS")
