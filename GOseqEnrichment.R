

## To find genes responding to the Treatment at each timepoint:

library(readr)
library(topGO)
library(dplyr)
library(goseq)
library(ggplot2)
library(tidyverse)

# for individual contrasts load: 
fit <- readRDS("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/Ind_fit.RDS")
my.contrasts <- readRDS("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/IndContrasts.RDS")

# for grouped contrasts load:
fit <- readRDS("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/Group_fit.RDS")
my.contrasts <- readRDS("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/GroupContrasts.RDS")
qlf_Control <- glmQLFTest(fit, contrast=my.contrasts[,"control"])
summary(decideTests(qlf_Control))

# Within the Eurasian genotype
qlf_E12 <- glmQLFTest(fit, contrast=my.contrasts[,"E_12"])
summary(decideTests(qlf_E12))
qlf_E24 <- glmQLFTest(fit, contrast=my.contrasts[,"E_24"])
summary(decideTests(qlf_E24))
qlf_E48 <- glmQLFTest(fit, contrast=my.contrasts[,"E_48"])
summary(decideTests(qlf_E48))
qlf_E96 <- glmQLFTest(fit, contrast=my.contrasts[,"E_96"])
summary(decideTests(qlf_E96))

# Within the hybrid genotype
qlf_H12 <- glmQLFTest(fit, contrast=my.contrasts[,"H_12"])
summary(decideTests(qlf_H12))
qlf_H24 <- glmQLFTest(fit, contrast=my.contrasts[,"H_24"])
summary(decideTests(qlf_H24))
qlf_H48 <- glmQLFTest(fit, contrast=my.contrasts[,"H_48"])
summary(decideTests(qlf_H48))
qlf_H96 <- glmQLFTest(fit, contrast=my.contrasts[,"H_96"])
summary(decideTests(qlf_H96))

# Between genotypes
qlf_GENO12 <- glmQLFTest(fit, contrast=my.contrasts[,"GENO_12"])
summary(decideTests(qlf_GENO12))
qlf_GENO24 <- glmQLFTest(fit, contrast=my.contrasts[,"GENO_24"])
summary(decideTests(qlf_GENO24))
qlf_GENO48 <- glmQLFTest(fit, contrast=my.contrasts[,"GENO_48"])
summary(decideTests(qlf_GENO48))
qlf_GENO96 <- glmQLFTest(fit, contrast=my.contrasts[,"GENO_96"])
summary(decideTests(qlf_GENO96))

## Build data files needed for enrichment analysis

txi <- readRDS("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/TxImport.RDS")

# Get Gene Ontology file in the right format for GOseq

go_data <- topGO::readMappings("/Users/Greg/Desktop/TrinotateTranscriptome/gene_ontology_proc.txt")
n_go <- lapply(go_data, length) %>% unlist
gene_go <- data.frame(gene_id=rep(names(go_data), n_go),
                      go_term=unlist(go_data)) %>%
  arrange(gene_id) %>%
  unique()
head(gene_go)

# Build average Length and Count DB 

avg_len <- as.data.frame(rowMeans(txi$length))
avg_count <- as.data.frame(rowMeans(txi$counts))
len_count_data <- merge(x = avg_len, y = avg_count, by = "row.names")

## GOseq for Control

# build 0/1 vector of all DEGs by merging with DB above to make sure everything is in the same order
Control_DE <- as.data.frame(abs(decideTests(qlf_Control)))
Control_DE <- merge(x = Control_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
Control_vector <- as.vector(Control_DE$`1*E_C -1*H_C`)
names(Control_vector) = Control_DE$Row.names

# Calculate bias using length
Control_pwf <- nullp(DEgenes = Control_vector, genome = "NULL", id = "Null", bias.data = Control_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin Hochberg 1995)
Control_GO.wall <- goseq(pwf = Control_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(Control_GO.wall)
Control_enriched.GO <- Control_GO.wall$category[p.adjust(Control_GO.wall$over_represented_pvalue, method = "BH")<.05]
Control_enriched.GO <- subset(Control_GO.wall, category %in% Control_enriched.GO)
head(Control_enriched.GO)

# Print top 10 enriched terms definitions
for(go in Control_enriched.GO$category[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n") 
}

Control_enriched.GO %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "num DEGs", fill = "P-value", title = "Control All DEGs Enrichment") +
  coord_flip()



## GOseq for E12

# build 0/1 vector of all DEGs by merging with DB above to make sure everything is in the same order
E12_DE <- as.data.frame(abs(decideTests(qlf_E12)))
E12_DE <- merge(x = E12_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
E12_vector <- as.vector(E12_DE$`-1*E_C_12 1*E_T_12`)
names(E12_vector) = E12_DE$Row.names

# Calculate bias using length
E12_pwf <- nullp(DEgenes = E12_vector, genome = "NULL", id = "Null", bias.data = E12_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin Hochberg 1995)
E12_GO.wall <- goseq(pwf = E12_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(E12_GO.wall)
E12_enriched.GO <- E12_GO.wall$category[p.adjust(E12_GO.wall$over_represented_pvalue, method = "BH")<.05]
E12_enriched.GO <- subset(E12_GO.wall, category %in% E12_enriched.GO)
head(E12_enriched.GO)

# Print top 10 enriched terms definitions
for(go in E12_enriched.GO$category[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n") 
}

E12_enriched.GO %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "num DEGs", fill = "P-value", title = "E12 All DEGs Enrichment") +
  coord_flip()


## GOseq for E24


# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
E24_DE <- as.data.frame(abs(decideTests(qlf_E24)))
E24_DE <- merge(x = E24_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
E24_vector <- as.vector(E24_DE$`-1*E_C_24 1*E_T_24`)
names(E24_vector) = E24_DE$Row.names 

# Calculate bias using length of gene
E24_pwf <- nullp(DEgenes = E24_vector, genome = "NULL", id = "Null", bias.data = E24_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin Hochberg 1995)
E24_GO.wall <- goseq(pwf = E24_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(E24_GO.wall)
E24_enriched.GO <- E24_GO.wall$category[p.adjust(E24_GO.wall$over_represented_pvalue, method = "BH")<.05]
E24_enriched.GO <- subset(E24_GO.wall, category %in% E24_enriched.GO)
head(E24_enriched.GO)

# Print top 10 enriched terms definitions
for(go in E24_enriched.GO$category[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n") 
}

E24_enriched.GO %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "num DE genes", fill = "P-value") +
  coord_flip()



## GOseq for E48



# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
E48_DE <- as.data.frame(abs(decideTests(qlf_E48)))
E48_DE <- merge(x = E48_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
E48_vector <- as.vector(E48_DE$`-1*E_C_48 1*E_T_48`)
names(E48_vector) = E48_DE$Row.names 

# Calculate bias using length of gene
E48_pwf <- nullp(DEgenes = E48_vector, genome = "NULL", id = "Null", bias.data = E48_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin Hochberg 1995)
E48_GO.wall <- goseq(pwf = E48_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(E48_GO.wall)
E48_enriched.GO <- E48_GO.wall$category[p.adjust(E48_GO.wall$over_represented_pvalue, method = "BH")<.05]
E48_enriched.GO <- subset(E48_GO.wall, category %in% E48_enriched.GO)
head(E48_enriched.GO)

# Print top 10 enriched terms definitions
for(go in E48_enriched.GO$category[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n") 
}

E48_enriched.GO %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "num DE genes", fill = "P-value") +
  coord_flip()



## GOseq for E96


# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
E96_DE <- as.data.frame(abs(decideTests(qlf_E96)))
E96_DE <- merge(x = E96_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
E96_vector <- as.vector(E96_DE$`-1*E_C_96 1*E_T_96`)
names(E96_vector) = E96_DE$Row.names 

# Calculate bias using length of gene
E96_pwf <- nullp(DEgenes = E96_vector, genome = "NULL", id = "Null", bias.data = E96_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin Hochberg 1995)
E96_GO.wall <- goseq(pwf = E96_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(E48_GO.wall)
E96_enriched.GO <- E96_GO.wall$category[p.adjust(E48_GO.wall$over_represented_pvalue, method = "BH")<.05]
E96_enriched.GO <- subset(E96_GO.wall, category %in% E96_enriched.GO)
head(E96_enriched.GO)

# Print top 10 enriched terms definitions
for(go in E96_enriched.GO$category[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n") 
}

E96_enriched.GO %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "num DE genes", fill = "P-value") +
  coord_flip()




## GOseq for H12



# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
H12_DE <- as.data.frame(abs(decideTests(qlf_H12)))
H12_DE <- merge(x = H12_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
H12_vector <- as.vector(H12_DE$`-1*H_C_12 1*H_T_12`)
names(H12_vector) = H12_DE$Row.names 

# Calculate bias using length
H12_pwf <- nullp(DEgenes = H12_vector, genome = "NULL", id = "Null", bias.data = H12_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin Hochberg 1995)
H12_GO.wall <- goseq(pwf = H12_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(H12_GO.wall)
H12_enriched.GO <- H12_GO.wall$category[p.adjust(H12_GO.wall$over_represented_pvalue, method = "BH")<.05]
H12_enriched.GO <- subset(H12_GO.wall, category %in% H12_enriched.GO)
head(H12_enriched.GO)

# Print top 10 enriched terms definitions
for(go in H12_enriched.GO$category[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n") 
}

H12_enriched.GO %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "num DE genes", fill = "P-value") +
  coord_flip()




## GOseq for H24



# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
H24_DE <- as.data.frame(abs(decideTests(qlf_H24)))
H24_DE <- merge(x = H24_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
H24_vector <- as.vector(H24_DE$`-1*H_C_24 1*H_T_24`)
names(H24_vector) = H24_DE$Row.names 

# Calculate bias using length of gene
H24_pwf <- nullp(DEgenes = H24_vector, genome = "NULL", id = "Null", bias.data = H24_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin Hochberg 1995)
H24_GO.wall <- goseq(pwf = H24_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(H24_GO.wall)
H24_enriched.GO <- H24_GO.wall$category[p.adjust(H24_GO.wall$over_represented_pvalue, method = "BH")<.05]
H24_enriched.GO <- subset(H24_GO.wall, category %in% H24_enriched.GO)
head(H24_enriched.GO)

# Print top 10 enriched terms definitions
for(go in H24_enriched.GO$category[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n") 
}

H24_enriched.GO %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "num DE genes", fill = "P-value") +
  coord_flip()



## GOseq for H48



# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
H48_DE <- as.data.frame(abs(decideTests(qlf_H48)))
H48_DE <- merge(x = H48_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
H48_vector <- as.vector(H48_DE$`-1*H_C_48 1*H_T_48`)
names(H48_vector) = H48_DE$Row.names 

# Calculate bias using length of gene
H48_pwf <- nullp(DEgenes = H48_vector, genome = "NULL", id = "Null", bias.data = H48_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin Hochberg 1995)
H48_GO.wall <- goseq(pwf = H48_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(H48_GO.wall)
H48_enriched.GO <- H48_GO.wall$category[p.adjust(H48_GO.wall$over_represented_pvalue, method = "BH")<.05]
H48_enriched.GO <- subset(H48_GO.wall, category %in% H48_enriched.GO)
head(H48_enriched.GO)

# Print top 10 enriched terms definitions
for(go in H48_enriched.GO$category[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n") 
}

H48_enriched.GO %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "num DE genes", fill = "P-value") +
  coord_flip()



## GOseq for H96



# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
H96_DE <- as.data.frame(abs(decideTests(qlf_H96)))
H96_DE <- merge(x = H96_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
H96_vector <- as.vector(H96_DE$`-1*H_C_96 1*H_T_96`)
names(H96_vector) = H96_DE$Row.names 

# Calculate bias using length of gene
H96_pwf <- nullp(DEgenes = H96_vector, genome = "NULL", id = "Null", bias.data = H96_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin Hochberg 1995)
H96_GO.wall <- goseq(pwf = H96_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(H96_GO.wall)
H96_enriched.GO <- H96_GO.wall$category[p.adjust(H96_GO.wall$over_represented_pvalue, method = "BH")<.05]
H96_enriched.GO <- subset(H96_GO.wall, category %in% H96_enriched.GO)
head(H96_enriched.GO)

# Print top 10 enriched terms definitions
for(go in H96_enriched.GO$category[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n") 
}

H96_enriched.GO %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "num DE genes", fill = "P-value") +
  coord_flip()


## GOseq for E96 vs H96



# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
G96_DE <- as.data.frame(abs(decideTests(qlf_GENO96)))
G96_DE <- merge(x = G96_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
G96_vector <- as.vector(G96_DE$`-1*E_C_96 1*E_T_96 1*H_C_96 -1*H_T_96`)
names(G96_vector) = G96_DE$Row.names 

# Calculate bias using length of gene
G96_pwf <- nullp(DEgenes = G96_vector, genome = "NULL", id = "Null", bias.data = G96_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin Hochberg 1995)
G96_GO.wall <- goseq(pwf = G96_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(G96_GO.wall)
G96_enriched.GO <- G96_GO.wall$category[p.adjust(G96_GO.wall$over_represented_pvalue, method = "BH")<.05]
G96_enriched.GO <- subset(G96_GO.wall, category %in% G96_enriched.GO)
head(G96_enriched.GO)

G96_under.GO <- G96_GO.wall$category[p.adjust(G96_GO.wall$under_represented_pvalue, method = "BH")<.05]
G96_under.GO <- subset(G96_GO.wall, category %in% G96_under.GO)
head(G96_under.GO)

# Print top 10 enriched terms definitions
for(go in G96_enriched.GO$category[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n") 
}

# Print top 10 under enriched terms definitions
for(go in G96_under.GO$category[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n") 
}

G96_enriched.GO %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "num DE genes", fill = "P-value") +
  theme(axis.text.x = element_text(size = 6)) +
  coord_flip()

G96_under.GO %>%
  arrange(desc(under_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat, fill = under_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "num DE genes", fill = "P-value") +
  coord_flip()



## GOseq for E48 vs H48


# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
G48_DE <- as.data.frame(abs(decideTests(qlf_GENO48)))
G48_DE <- merge(x = G48_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
G48_vector <- as.vector(G48_DE$`-1*E_C_48 1*E_T_48 1*H_C_48 -1*H_T_48`)
names(G48_vector) = G48_DE$Row.names 

# Calculate bias using length of gene
G48_pwf <- nullp(DEgenes = G48_vector, genome = "NULL", id = "Null", bias.data = G48_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin Hochberg 1995)
G48_GO.wall <- goseq(pwf = G48_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(G48_GO.wall)
G48_enriched.GO <- G48_GO.wall$category[p.adjust(G48_GO.wall$over_represented_pvalue, method = "BH")<.05]
G48_enriched.GO <- subset(G48_GO.wall, category %in% G48_enriched.GO)
head(G48_enriched.GO)

G48_under.GO <- G48_GO.wall$category[p.adjust(G48_GO.wall$under_represented_pvalue, method = "BH")<.05]
G48_under.GO <- subset(G48_GO.wall, category %in% G48_under.GO)
head(G48_under.GO)

# Print top 10 enriched terms definitions
for(go in G48_enriched.GO$category[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n") 
}

# Print top 10 under enriched terms definitions
for(go in G48_under.GO$category[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n") 
}

G48_enriched.GO %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "num DE genes", fill = "P-value") +
  coord_flip()

G48_under.GO %>%
  arrange(desc(under_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat, fill = under_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "num DE genes", fill = "P-value") +
  coord_flip()



## Venn Diagram plotting of numbers of shared DEGs between contrasts



library(VennDetail)
E12ven.vec <- E12_DE$Row.names[E12_DE$`-1*E_C_12 1*E_T_12`==1]
E24ven.vec <- E24_DE$Row.names[E24_DE$`-1*E_C_24 1*E_T_24`==1]
E48ven.vec <- E48_DE$Row.names[E48_DE$`-1*E_C_48 1*E_T_48`==1]
E96ven.vec <- E96_DE$Row.names[E96_DE$`-1*E_C_96 1*E_T_96`==1]

H12ven.vec <- H12_DE$Row.names[H12_DE$`-1*H_C_12 1*H_T_12`==1]
H24ven.vec <- H24_DE$Row.names[H24_DE$`-1*H_C_24 1*H_T_24`==1]
H48ven.vec <- H48_DE$Row.names[H48_DE$`-1*H_C_48 1*H_T_48`==1]
H96ven.vec <- H96_DE$Row.names[H96_DE$`-1*H_C_96 1*H_T_96`==1]

E_ven <- venndetail(list(E12 = E12ven.vec, E24 = E24ven.vec, E48 = E48ven.vec, E96 = E96ven.vec))
H_ven <- venndetail(list(H12 = H12ven.vec, H24 = H24ven.vec, H48 = H48ven.vec, H96 = H96ven.vec))
Early_ven <- venndetail(list(E12 = E12ven.vec, E24 = E24ven.vec, H12 = H12ven.vec, H24 = H24ven.vec))
ven12 <- venndetail(list(E12 = E12ven.vec, H12 = H12ven.vec))
ven24 <- venndetail(list(E24 = E24ven.vec, H24 = H24ven.vec))
ven48 <- venndetail(list(E48 = E48ven.vec, H48 = H48ven.vec))
ven96 <- venndetail(list(E96 = E96ven.vec, H96 = H96ven.vec))
plot(ven12)
graphics.off()
plot(ven24)
graphics.off()
plot(ven48)
graphics.off()
plot(ven96)
graphics.off()
plot(E_ven)
graphics.off()
plot(H_ven)
graphics.off()
plot(Early_ven)
graphics.off()
