

## To find genes responding to the Treatment at each timepoint:

library(edgeR)
library(readr)
library(topGO)
library(dplyr)
library(goseq)
library(ggplot2)
library(tidyverse)
library(data.table)
library(VennDetail)

# for individual contrasts load: 
#fit <- readRDS("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/Ind_fit.RDS")
#my.contrasts <- readRDS("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/IndContrasts.RDS")

# for grouped contrasts load:
fit <- readRDS("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/Group_fit.RDS")
my.contrasts <- readRDS("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/GroupContrasts.RDS")

# Perform contrasts between groups for DE and write the DE table to a file.

# Between the control groups (How is expression constitutively different between genotypes?)
qlf_Control <- glmQLFTest(fit, contrast=my.contrasts[,"control"])
write.csv(qlf_Control$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_control.csv")


# Within the Eurasian genotype over time (How is expression affected by treatment at each of the time points after treatment?)
qlf_E12 <- glmQLFTest(fit, contrast=my.contrasts[,"E_12"])
write.csv(qlf_E12$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_E12.csv")

qlf_E24 <- glmQLFTest(fit, contrast=my.contrasts[,"E_24"])
write.csv(qlf_E24$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_E24.csv")

qlf_E48 <- glmQLFTest(fit, contrast=my.contrasts[,"E_48"])
write.csv(qlf_E48$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_E48.csv")

qlf_E96 <- glmQLFTest(fit, contrast=my.contrasts[,"E_96"])
write.csv(qlf_E96$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_E96.csv")


# Within the hybrid genotype over time (How is expression affected by treatment at each of the time points after treatment?)
qlf_H12 <- glmQLFTest(fit, contrast=my.contrasts[,"H_12"])
write.csv(qlf_H12$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_H12.csv")

qlf_H24 <- glmQLFTest(fit, contrast=my.contrasts[,"H_24"])
write.csv(qlf_H24$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_H24.csv")

qlf_H48 <- glmQLFTest(fit, contrast=my.contrasts[,"H_48"])
write.csv(qlf_H48$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_H48.csv")

qlf_H96 <- glmQLFTest(fit, contrast=my.contrasts[,"H_96"])
write.csv(qlf_H96$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_H96.csv")

qlfSummary <- data.table(
  E12 = summary(decideTests(qlf_E12)), H12 = summary(decideTests(qlf_H12)), 
  E24 = summary(decideTests(qlf_E24)), H24 = summary(decideTests(qlf_H24)),
  E48 = summary(decideTests(qlf_E48)), H48 = summary(decideTests(qlf_H48)),
  E96 = summary(decideTests(qlf_E96)), H96 = summary(decideTests(qlf_H96)))

write_csv(qlfSummary, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/updown.csv")

# Interactions between the Eurasian and hybrid (How is expression affected by treatment at each of the time points after treatment differently between the 2 genotypes?)
qlf_HAT12 <- glmQLFTest(fit, contrast=my.contrasts[,"HAT_12"])
write.csv(qlf_HAT12$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_HAT12.csv")

qlf_HAT24 <- glmQLFTest(fit, contrast=my.contrasts[,"HAT_24"])
write.csv(qlf_HAT24$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_HAT24.csv")

qlf_HAT48 <- glmQLFTest(fit, contrast=my.contrasts[,"HAT_48"])
write.csv(qlf_HAT48$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_HAT48.csv")

qlf_HAT96 <- glmQLFTest(fit, contrast=my.contrasts[,"HAT_96"])
write.csv(qlf_HAT96$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_HAT96.csv")

# Between treatment and control within genotypes (How does treatment affect expression as a whole?)
qlf_E_treat <- glmQLFTest(fit, contrast=my.contrasts[,"E_treat"])
write.csv(qlf_E_treat$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_E_treat.csv")

qlf_H_treat <- glmQLFTest(fit, contrast=my.contrasts[,"H_treat"])
write.csv(qlf_H_treat$table, row.names = TRUE, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/qlf_H_treat.csv")

qlf_HAT_Summary <- data.table(
  HAT12 = summary(decideTests(qlf_HAT12)), 
  HAT24 = summary(decideTests(qlf_HAT24)),
  HAT48 = summary(decideTests(qlf_HAT48)),
  HAT96 = summary(decideTests(qlf_HAT96)))

write_csv(qlf_HAT_Summary, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/updown_interaction.csv")

## Build data files needed for enrichment analysis

txi <- readRDS("/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/TxImport.RDS")

# Get Gene Ontology file in the right format for GOseq

go_data <- topGO::readMappings("/Users/Greg/Documents/GitHub/WMF_DeNovoTranscriptome/results/TrinotateTranscriptome/gene_ontology_proc.txt")
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



## GOseq for Controls

# build 0/1 vector of all DEGs by merging with DB above to make sure everything is in the same order
Control_DE <- as.data.frame(abs(decideTests(qlf_Control)))
Control_DE <- merge(x = Control_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
Control_vector <- as.vector(Control_DE[,2])
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
#for(go in Control_enriched.GO$category[1:10]){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n") }

write_csv(Control_enriched.GO, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/Control.enrichment.csv")

Control_enriched.GO %>%
  filter(ontology == 'CC') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "Control CC") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("Control.CC.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

Control_enriched.GO %>%
  filter(ontology == 'BP') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "Control BP") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("Control.BP.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

Control_enriched.GO %>%
  filter(ontology == 'MF') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "Control MF") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("Control.MF.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)




## GOseq for E12

# build 0/1 vector of all DEGs by merging with DB above to make sure everything is in the same order
E12_DE <- as.data.frame(abs(decideTests(qlf_E12)))
E12_DE <- merge(x = E12_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
E12_vector <- as.vector(E12_DE[,2])
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
#for(go in E12_enriched.GO$category[1:10]){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n") }

write_csv(E12_enriched.GO, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/E12.enrichment.csv")


E12_enriched.GO %>%
  filter(ontology == 'CC') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = ontology , fill = over_represented_pvalue, size=numDEInCat/numInCat)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E12 CC") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()


E12_enriched.GO %>%
  filter(ontology == 'BP') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E12 BP") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("E12.BP.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

E12_enriched.GO %>%
  filter(ontology == 'MF') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E12 MF") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("E12.MF.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)



## GOseq for E24


# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
E24_DE <- as.data.frame(abs(decideTests(qlf_E24)))
E24_DE <- merge(x = E24_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
E24_vector <- as.vector(E24_DE[,2])
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
#for(go in E24_enriched.GO$category[1:10]){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n") }


write_csv(E24_enriched.GO, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/E24.enrichment.csv")


E24_enriched.GO %>%
  filter(ontology == 'CC') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E24 CC") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("E24.CC.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

E24_enriched.GO %>%
  filter(ontology == 'BP') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E24 BP") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("E24.BP.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

E24_enriched.GO %>%
  filter(ontology == 'MF') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E24 MF") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("E24.MF.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)



## GOseq for E48



# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
E48_DE <- as.data.frame(abs(decideTests(qlf_E48)))
E48_DE <- merge(x = E48_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
E48_vector <- as.vector(E48_DE[,2])
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
#for(go in E48_enriched.GO$category[1:10]){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n") }

write_csv(E48_enriched.GO, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/E48.enrichment.csv")


E48_enriched.GO %>%
  filter(ontology == 'CC') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E48 CC") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("E48.CC.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

E48_enriched.GO %>%
  filter(ontology == 'BP') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E48 BP") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("E48.BP.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

E48_enriched.GO %>%
  filter(ontology == 'MF') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E48 MF") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("E48.MF.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)




## GOseq for E96


# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
E96_DE <- as.data.frame(abs(decideTests(qlf_E96)))
E96_DE <- merge(x = E96_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
E96_vector <- as.vector(E96_DE[,2])
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
#for(go in E96_enriched.GO$category[1:10]){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n") }

write_csv(E96_enriched.GO, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/E96.enrichment.csv")


E96_enriched.GO %>%
  filter(ontology == 'CC') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E96 CC") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("E96.CC.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

E96_enriched.GO %>%
  filter(ontology == 'BP') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E96 BP") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("E96.BP.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

E96_enriched.GO %>%
  filter(ontology == 'MF') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E96 MF") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("E96.MF.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)





## GOseq for H12


# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
H12_DE <- as.data.frame(abs(decideTests(qlf_H12)))
H12_DE <- merge(x = H12_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
H12_vector <- as.vector(H12_DE[,2])
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
#for(go in H12_enriched.GO$category[1:10]){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n") }

write_csv(H12_enriched.GO, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/H12.enrichment.csv")


H12_enriched.GO %>%
  filter(ontology == 'CC') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H12 CC") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H12.CC.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

H12_enriched.GO %>%
  filter(ontology == 'BP') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H12 BP") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H12.BP.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

H12_enriched.GO %>%
  filter(ontology == 'MF') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H12 MF") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H12.MF.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)



## GOseq for H24



# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
H24_DE <- as.data.frame(abs(decideTests(qlf_H24)))
H24_DE <- merge(x = H24_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
H24_vector <- as.vector(H24_DE[,2])
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
#for(go in H24_enriched.GO$category[1:10]){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n") }

write_csv(H24_enriched.GO, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/H24.enrichment.csv")


H24_enriched.GO %>%
  filter(ontology == 'CC') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H24 CC") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H24.CC.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

H24_enriched.GO %>%
  filter(ontology == 'BP') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H24 BP") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H24.BP.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

H24_enriched.GO %>%
  filter(ontology == 'MF') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H24 MF") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H24.MF.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)



## GOseq for H48



# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
H48_DE <- as.data.frame(abs(decideTests(qlf_H48)))
H48_DE <- merge(x = H48_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
H48_vector <- as.vector(H48_DE[,2])
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
#for(go in H48_enriched.GO$category[1:10]){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n") }

write_csv(H48_enriched.GO, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/H48.enrichment.csv")


H48_enriched.GO %>%
  filter(ontology == 'CC') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H48 CC") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H48.CC.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

H48_enriched.GO %>%
  filter(ontology == 'BP') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H48 BP") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H48.BP.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

H48_enriched.GO %>%
  filter(ontology == 'MF') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H48 MF") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H48.MF.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)



## GOseq for H96


# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
H96_DE <- as.data.frame(abs(decideTests(qlf_H96)))
H96_DE <- merge(x = H96_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
H96_vector <- as.vector(H96_DE[,2])
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
#for(go in H96_enriched.GO$category[1:10]){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n") }

write_csv(H96_enriched.GO, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/H96.enrichment.csv")


H96_enriched.GO %>%
  filter(ontology == 'CC') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H96 CC") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H96.CC.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

H96_enriched.GO %>%
  filter(ontology == 'BP') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H96 BP") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H96.BP.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

H96_enriched.GO %>%
  filter(ontology == 'MF') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H96 MF") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H96.MF.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)


## GOseq for HAT12 **NO INTERACTION GENES AT 12HAT**


## GOseq for HAT24 **Only 3 genes in up and 3 in down = Not enough to calculate enrichment**


## GOseq for HAT48

# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
HAT48_DE <- as.data.frame(abs(decideTests(qlf_HAT48)))
HAT48_DE <- merge(x = HAT48_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
HAT48_vector <- as.vector(HAT48_DE[,2])
names(HAT48_vector) = HAT48_DE$Row.names 

# Calculate bias using length of gene
HAT48_pwf <- nullp(DEgenes = HAT48_vector, genome = "NULL", id = "Null", bias.data = HAT48_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin HATochberg 1995)
HAT48_GO.wall <- goseq(pwf = HAT48_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(HAT48_GO.wall)
HAT48_enriched.GO <- HAT48_GO.wall$category[p.adjust(HAT48_GO.wall$over_represented_pvalue, method = "BH")<.05]
HAT48_enriched.GO <- subset(HAT48_GO.wall, category %in% HAT48_enriched.GO)
head(HAT48_enriched.GO)

# Print top 10 enriched terms definitions
#for(go in HAT48_enriched.GO$category[1:10]){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n") }

write_csv(HAT48_enriched.GO, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/HAT48.enrichment.csv")


HAT48_enriched.GO %>%
  filter(ontology == 'CC') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "HAT48 CC") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("HAT48.CC.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

HAT48_enriched.GO %>%
  filter(ontology == 'BP') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "HAT48 BP") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("HAT48.BP.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

HAT48_enriched.GO %>%
  filter(ontology == 'MF') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "HAT48 MF") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("HAT48.MF.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)



## GOseq for HAT96


# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
HAT96_DE <- as.data.frame(abs(decideTests(qlf_HAT96)))
HAT96_DE <- merge(x = HAT96_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
HAT96_vector <- as.vector(HAT96_DE[,2])
names(HAT96_vector) = HAT96_DE$Row.names 

# Calculate bias using length of gene
HAT96_pwf <- nullp(DEgenes = HAT96_vector, genome = "NULL", id = "Null", bias.data = HAT96_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin HATochberg 1995)
HAT96_GO.wall <- goseq(pwf = HAT96_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(HAT96_GO.wall)
HAT96_enriched.GO <- HAT96_GO.wall$category[p.adjust(HAT96_GO.wall$over_represented_pvalue, method = "BH")<.05]
HAT96_enriched.GO <- subset(HAT96_GO.wall, category %in% HAT96_enriched.GO)
head(HAT96_enriched.GO)

# Print top 10 enriched terms definitions
#for(go in HAT96_enriched.GO$category[1:10]){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n") }

write_csv(HAT96_enriched.GO, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/HAT96.enrichment.csv")


HAT96_enriched.GO %>%
  filter(ontology == 'CC') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "HAT96 CC") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("HAT96.CC.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

HAT96_enriched.GO %>%
  filter(ontology == 'BP') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "HAT96 BP") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("HAT96.BP.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

HAT96_enriched.GO %>%
  filter(ontology == 'MF') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "HAT96 MF") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("HAT96.MF.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

## GOseq for E_treat


# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
E_treat_DE <- as.data.frame(abs(decideTests(qlf_E_treat)))
E_treat_DE <- merge(x = E_treat_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
E_treat_vector <- as.vector(E_treat_DE[,2])
names(E_treat_vector) = E_treat_DE$Row.names 

# Calculate bias using length of gene
E_treat_pwf <- nullp(DEgenes = E_treat_vector, genome = "NULL", id = "Null", bias.data = E_treat_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin Hochberg 1995)
E_treat_GO.wall <- goseq(pwf = E_treat_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(E_treat_GO.wall)
E_treat_enriched.GO <- E_treat_GO.wall$category[p.adjust(E_treat_GO.wall$over_represented_pvalue, method = "BH")<.05]
E_treat_enriched.GO <- subset(E_treat_GO.wall, category %in% E_treat_enriched.GO)
head(E_treat_enriched.GO)

# Print top 10 enriched terms definitions
#for(go in E_treat_enriched.GO$category[1:10]){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n") }

write_csv(E_treat_enriched.GO, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/E_treat.enrichment.csv")


E_treat_enriched.GO %>%
  filter(ontology == 'CC') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E_treat CC") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("E_treat.CC.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

E_treat_enriched.GO %>%
  filter(ontology == 'BP') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E_treat BP") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("E_treat.BP.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

E_treat_enriched.GO %>%
  filter(ontology == 'MF') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "E_treat MF") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("E_treat.MF.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)



## GOseq for H_treat


# build 0/1 vector of DE by merging with DB above to make sure everything is in the same order
H_treat_DE <- as.data.frame(abs(decideTests(qlf_H_treat)))
H_treat_DE <- merge(x = H_treat_DE, y = len_count_data, by.x = "row.names", by.y = "Row.names")
H_treat_vector <- as.vector(H_treat_DE[,2])
names(H_treat_vector) = H_treat_DE$Row.names 

# Calculate bias using length of gene
H_treat_pwf <- nullp(DEgenes = H_treat_vector, genome = "NULL", id = "Null", bias.data = H_treat_DE$`rowMeans(txi$length)`)

# Find over and under enrichment Go terms with P-value correction (0.05 FDR with Benjamin Hochberg 1995)
H_treat_GO.wall <- goseq(pwf = H_treat_pwf, genome = "NULL", id = "NULL", gene2cat = gene_go)
head(H_treat_GO.wall)
H_treat_enriched.GO <- H_treat_GO.wall$category[p.adjust(H_treat_GO.wall$over_represented_pvalue, method = "BH")<.05]
H_treat_enriched.GO <- subset(H_treat_GO.wall, category %in% H_treat_enriched.GO)
head(H_treat_enriched.GO)

# Print top 10 enriched terms definitions
#for(go in H_treat_enriched.GO$category[1:10]){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n") }

write_csv(H_treat_enriched.GO, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/H_treat.enrichment.csv")


H_treat_enriched.GO %>%
  filter(ontology == 'CC') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H_treat CC") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H_treat.CC.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

H_treat_enriched.GO %>%
  filter(ontology == 'BP') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H_treat BP") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H_treat.BP.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)

H_treat_enriched.GO %>%
  filter(ontology == 'MF') %>%
  arrange(desc(over_represented_pvalue)) %>%
  ggplot(aes(x = fct_inorder(term), y = numDEInCat/numInCat, fill = over_represented_pvalue)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
  geom_bar(stat = "identity") + 
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  labs(x = "", y = "Propotion of catalog DE", fill = "P-value", title = "H_treat MF") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

ggsave("H_treat.MF.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 320)


### Grouped GO enrichment plots by timepoint

E12_enriched.GO <- E12_enriched.GO %>%
  mutate(genotype = "E12")
H12_enriched.GO <- H12_enriched.GO %>%
  mutate(genotype = "H12")

GO.12 <- bind_rows(E12_enriched.GO, H12_enriched.GO)

GO.12.compare <- E12_enriched.GO %>%
  full_join(H12_enriched.GO, by = 'category', copy = TRUE, suffix = c("",".H12"))

write_csv(GO.12.compare, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/GO.12.enrichment.csv")

GO.12 %>%
    arrange(desc(over_represented_pvalue)) %>%
    na.omit() %>%
    ggplot(aes(x = fct_inorder(term), y = genotype)) + 
    geom_point(aes(color = over_represented_pvalue, size = numDEInCat/numInCat, group = ontology), alpha = 1) + 
    scale_color_gradient(low = "blue", high = "red", na.value = NA) +
    scale_x_discrete(label = function(x) stringr::str_trunc(x, 40))+
    labs(x = "", y = "", color = "P-value", size = "ratio of DE to catalog") +
    theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 0)) +
    coord_flip() +
    facet_grid(rows = vars(ontology), scale = "free", space = "free")


ggsave("GO.12.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 600,
  width = 9,
  height = 16,
  units = "in")

GO.12 %>%
    filter(ontology == 'MF') %>%
    arrange(desc(over_represented_pvalue)) %>%
    na.omit() %>%
    ggplot(aes(x = fct_inorder(term), y = genotype)) + 
    geom_point(aes(color = over_represented_pvalue, size = numDEInCat/numInCat, group = ontology), alpha = 1) + 
    scale_color_gradient(low = "blue", high = "red", na.value = NA) +
    scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
    labs(x = "", y = "", color = "P-value", size = "ratio of DE to catalog") +
    theme(axis.text.y = element_text(size = 8)) +
    coord_flip()

ggsave("GO.12.MF.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 600,
  width = 9,
  height = 16,
  units = "in")

GO.12 %>%
    filter(ontology == 'BP') %>%
    arrange(desc(over_represented_pvalue)) %>%
    na.omit() %>%
    ggplot(aes(x = fct_inorder(term), y = genotype)) + 
    geom_point(aes(color = over_represented_pvalue, size = numDEInCat/numInCat, group = ontology), alpha = 1) + 
    scale_color_gradient(low = "blue", high = "red", na.value = NA) +
    scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
    labs(x = "", y = "", color = "P-value", size = "ratio of DE to catalog") +
    theme(axis.text.y = element_text(size = 8)) +
    coord_flip()

ggsave("GO.12.BP.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 600,
  width = 9,
  height = 16,
  units = "in")

E24_enriched.GO <- E24_enriched.GO %>%
  mutate(genotype = "E24")
H24_enriched.GO <- H24_enriched.GO %>%
  mutate(genotype = "H24")

GO.24 <- bind_rows(E24_enriched.GO, H24_enriched.GO)


GO.24.compare <- E24_enriched.GO %>%
  full_join(H24_enriched.GO, by = 'category', copy = TRUE, suffix = c("",".H24"))

write_csv(GO.24.compare, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/GO.24.enrichment.csv")

GO.24 %>%
    arrange(desc(over_represented_pvalue)) %>%
    na.omit() %>%
    ggplot(aes(x = fct_inorder(term), y = genotype)) + 
    geom_point(aes(color = over_represented_pvalue, size = numDEInCat/numInCat, group = ontology), alpha = 1) + 
    scale_color_gradient(low = "blue", high = "red", na.value = NA) +
    scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
    labs(x = "", y = "", color = "P-value", size = "ratio of DE to catalog") +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip() +
    facet_grid(rows = vars(ontology), scale = "free", space = "free")

ggsave("GO.24.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 600,
  width = 9,
  height = 16,
  units = "in")


E48_enriched.GO <- E48_enriched.GO %>%
  mutate(genotype = "E48")
H48_enriched.GO <- H48_enriched.GO %>%
  mutate(genotype = "H48")

GO.48 <- bind_rows(E48_enriched.GO, H48_enriched.GO)


GO.48.compare <- E48_enriched.GO %>%
  full_join(H48_enriched.GO, by = 'category', copy = TRUE, suffix = c("",".H48"))

write_csv(GO.48.compare, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/GO.48.enrichment.csv")

GO.48 %>%
    arrange(desc(over_represented_pvalue)) %>%
    na.omit() %>%
    ggplot(aes(x = fct_inorder(term), y = genotype)) + 
    geom_point(aes(color = over_represented_pvalue, size = numDEInCat/numInCat, group = ontology), alpha = 1) + 
    scale_color_gradient(low = "blue", high = "red", na.value = NA) +
    scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
    labs(x = "", y = "", color = "P-value", size = "ratio of DE to catalog") +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip() +
    facet_grid(rows = vars(ontology), scale = "free", space = "free")

ggsave("GO.48.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 600,
  width = 9,
  height = 19,
  units = "in")

E96_enriched.GO <- E96_enriched.GO %>%
  mutate(genotype = "E96")
H96_enriched.GO <- H96_enriched.GO %>%
  mutate(genotype = "H96")

GO.96 <- bind_rows(E96_enriched.GO, H96_enriched.GO)



GO.96.compare <- E96_enriched.GO %>%
  full_join(H96_enriched.GO, by = 'category', copy = TRUE, suffix = c("",".H96"))

write_csv(GO.96.compare, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/GO.96.enrichment.csv")

GO.96 %>%
    arrange(desc(over_represented_pvalue)) %>%
    na.omit() %>%
    ggplot(aes(x = fct_inorder(term), y = genotype)) + 
    geom_point(aes(color = over_represented_pvalue, size = numDEInCat/numInCat, group = ontology), alpha = 1) + 
    scale_color_gradient(low = "blue", high = "red", na.value = NA) +
    scale_x_discrete(label = function(x) stringr::str_trunc(x, 40)) +
    labs(x = "", y = "", color = "P-value", size = "ratio of DE to catalog") +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip() +
    facet_grid(rows = vars(ontology), scale = "free", space = "free")

ggsave("GO.96.pdf",
  plot = last_plot(),
  device = NULL,
  path = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/",
  scale = 1,
  dpi = 600,
  width = 9,
  height = 19,
  units = "in")


## Venn Diagram plotting of numbers of shared DEGs between contrasts


E12ven.vec <- E12_DE$Row.names[E12_DE[,2]==1]
E24ven.vec <- E24_DE$Row.names[E24_DE[,2]==1]
E48ven.vec <- E48_DE$Row.names[E48_DE[,2]==1]
E96ven.vec <- E96_DE$Row.names[E96_DE[,2]==1]

H12ven.vec <- H12_DE$Row.names[H12_DE[,2]==1]
H24ven.vec <- H24_DE$Row.names[H24_DE[,2]==1]
H48ven.vec <- H48_DE$Row.names[H48_DE[,2]==1]
H96ven.vec <- H96_DE$Row.names[H96_DE[,2]==1]

E_treatven.vec <- E_treat_DE$Row.names[E_treat_DE[,2]==1]
H_treatven.vec <- H_treat_DE$Row.names[H_treat_DE[,2]==1]

Treat_ven<- venndetail(list(E = E_treatven.vec, H = H_treatven.vec))
E_ven <- venndetail(list(E12 = E12ven.vec, E24 = E24ven.vec, E48 = E48ven.vec, E96 = E96ven.vec))
H_ven <- venndetail(list(H12 = H12ven.vec, H24 = H24ven.vec, H48 = H48ven.vec, H96 = H96ven.vec))
Early_ven <- venndetail(list(E12 = E12ven.vec, E24 = E24ven.vec, H12 = H12ven.vec, H24 = H24ven.vec))
ven12 <- venndetail(list(E12 = E12ven.vec, H12 = H12ven.vec))
ven24 <- venndetail(list(E24 = E24ven.vec, H24 = H24ven.vec))
ven48 <- venndetail(list(E48 = E48ven.vec, H48 = H48ven.vec))
ven96 <- venndetail(list(E96 = E96ven.vec, H96 = H96ven.vec))


pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/E_Ven.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
plot(E_ven)
dev.off()

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/H_Ven.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
plot(H_ven)
dev.off()

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/Ven12.pdf",   # The directory you want to save the file in
    width = 24, # The width of the plot in inches
    height = 24) # The height of the plot in inches
plot(ven12, col = c("white", "Gray"), cex = 0.5)
dev.off()

shared12 <- ven12[0:sum(ven12$Subset=="Shared"),2]
write_csv(as.data.frame(shared12), "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_12.csv")

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/Ven24.pdf",   # The directory you want to save the file in
    width = 24, # The width of the plot in inches
    height = 24) # The height of the plot in inches
plot(ven24)
dev.off()

shared24 <- ven24[0:sum(ven24$Subset=="Shared"),2]
write_csv(as.data.frame(shared24), "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_24.csv")

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/Ven48.pdf",   # The directory you want to save the file in
    width = 24, # The width of the plot in inches
    height = 24) # The height of the plot in inches
plot(ven48)
dev.off()

shared48 <- ven48[0:sum(ven48$Subset=="Shared"),2]
write_csv(as.data.frame(shared48), "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_48.csv")

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/Ven96.pdf",   # The directory you want to save the file in
    width = 24, # The width of the plot in inches
    height = 24) # The height of the plot in inches
plot(ven96)
dev.off()

shared96 <- ven96[0:sum(ven96$Subset=="Shared"),2]
write_csv(as.data.frame(shared96), "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_96.csv")

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/VenEWM.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
plot(E_ven)
dev.off()

sharedE <- E_ven[0:sum(E_ven$Subset=="Shared"),2]
write_csv(as.data.frame(sharedE), "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_E.csv")

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/VenHWM.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
plot(H_ven)
dev.off()

sharedH <- H_ven[0:sum(H_ven$Subset=="Shared"),2]
write_csv(as.data.frame(sharedH), "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_H.csv")

pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/EarlyVen.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
plot(Early_ven)
dev.off()

sharedEarly <- Early_ven[0:sum(Early_ven$Subset=="Shared"),2]
write_csv(as.data.frame(sharedEarly), "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_Early.csv")


pdf(file = "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/TreatVen.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
plot(Treat_ven, main = "Treatment across all time")
dev.off()

sharedTreat <- Treat_ven[0:sum(Treat_ven$Subset=="Shared"),2]
write_csv(as.data.frame(sharedTreat), "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/shared_Treat.csv")



## Export lists of significant DEGs for future analyses
Control.sig <- Control_pwf %>% filter(DEgenes != 0)
write.csv(Control.sig, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/Control_sig.csv")

E12.sig <- E12_pwf %>% filter(DEgenes != 0)
write.csv(E12.sig, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/E12_sig.csv")
E24.sig <- E24_pwf %>% filter(DEgenes != 0)
write.csv(E24.sig, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/E24_sig.csv")
E48.sig <- E48_pwf %>% filter(DEgenes != 0)
write.csv(E48.sig, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/E48_sig.csv")
E96.sig <- E96_pwf %>% filter(DEgenes != 0)
write.csv(E96.sig, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/E96_sig.csv")

H12.sig <- H12_pwf %>% filter(DEgenes != 0)
write.csv(H12.sig, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/H12_sig.csv")
H24.sig <- H24_pwf %>% filter(DEgenes != 0)
write.csv(H24.sig, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/H24_sig.csv")
H48.sig <- H48_pwf %>% filter(DEgenes != 0)
write.csv(H48.sig, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/H48_sig.csv")
H96.sig <- H96_pwf %>% filter(DEgenes != 0)
write.csv(H96.sig, "/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/H96_sig.csv")


# Redraw the Venn Diagrams made in VennDetail with modified text and color
ven.12 <- draw.pairwise.venn(area1 = 322, area2 = 225, cross.area = 141, category = c("EWM","hybrid"), fill = c("white","grey"), ext.text = FALSE, cex = 3, cat.pos = c(335, 25), cat.cex = 2)
grid.newpage()
ven.24 <- draw.pairwise.venn(area1 = 2655, area2 = 570, cross.area = 449, category = c("EWM","hybrid"), fill = c("white","grey"), ext.text = FALSE, cex = 3, cat.pos = c(335, 25), cat.cex = 2)
grid.newpage()
ven.48 <- draw.pairwise.venn(area1 = 13146, area2 = 2455, cross.area = 1944, category = c("EWM","hybrid"), fill = c("white","grey"), ext.text = FALSE, cex = 3, cat.pos = c(335, 25), cat.cex = 2)
grid.newpage()
ven.96 <- draw.pairwise.venn(area1 = 7948, area2 = 2521, cross.area = 2054, category = c("EWM","hybrid"), fill = c("white","grey"), ext.text = FALSE, cex = 3, cat.pos = c(335, 25), cat.cex = 2)

# Draw stacked bar graph of the total number of up- and downregulated genes
ggplot(DEcounts_long, aes(fill = Direction, y = DEGs, x = Genotype)) + 
    geom_bar(position = "stack", stat = "identity") + 
    theme_minimal() + 
    xlab("") +
    facet_wrap(~HAT, strip.position = "bottom", scales = "free_x", ncol = 4) +
    theme(text = element_text(size = 20)) +
    theme(axis.text.x = element_text(size = 8))    







