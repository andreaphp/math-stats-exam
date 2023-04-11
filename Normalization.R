#Description and effect of normalization techniques
library(tidyverse)
library(ggfortify)
library(preprocessCore)
library(edgeR)

#DATASET CompCODEr
browseVignettes("compcodeR")


#musa article
#   In MuSA tool, radiomic and genomic features can be normalized accordingly to five most common nor-
#   malization methods. Min–Max normalization, Z-score normalization, log2 normalization, upper quartile and
# whitening methods were included in this step. Min–max normalization is one of the most common method
# to normalize data. For every feature, the minimum value of a feature gets transformed into a 0, the maximum
# value gets transformed into a 1, and every other value is rescaled to lie within a range of 0 to 1 ^36
# .Regarding stand-
#   ardized z-score normalization, each feature was normalized as z = (x-−
#                                                                      x)/s, where x, −
# x and s are the feature, the
# meanand the standard deviation respectively ^37 . Features can also be standardized as Log-transformation (base
#                                                                                                           2), a constant value a = b – min(x) where b is 1 and x is the feature, was added to the data for handling negative
# values ^38 . The upper quartile normalization divides each read count by the 75th percentile of the read counts in
# its sample ^39 ; lastly, whitening normalization technique from the principle component analysis (PCA), is based
# on a linear transformation that converts a vector of random variables with a known covariance matrix into a
# set of new variables whose covariance is the identity matrix, meaning that they are uncorrelated and each have
# variance equal to one ^40 .*/

#****************************************************************************
# Problem description ---------------------------------------------------------
#Lately, normalization techniques are of great interest in data preprocessing especially for machine learning algorithms ^35
#In fact, overlooking feature normalization step may lead to individual features being over or underrepresented and
#eventually introduce bias into developed models and statistical analysis ^34

#****************************************************************************
# Synthetic Dataset Generation (Vignette CompCodeR)
B_625_625 <- generateSyntheticData(dataset = "B_625_625", n.vars = 12500, 
                                   samples.per.cond = 5, n.diffexp = 1250, 
                                   repl.id = 1, seqdepth = 1e7, 
                                   fraction.upregulated = 0.5, 
                                   between.group.diffdisp = FALSE, 
                                   filter.threshold.total = 1, 
                                   filter.threshold.mediancpm = 0, 
                                   fraction.non.overdispersed = 0, 
                                   output.file = "B_625_625_5spc_repl1.rds")


summarizeSyntheticDataSet(data.set = "B_625_625_5spc_repl1.rds", 
                          output.filename = "B_625_625_5spc_repl1_datacheck.html")

#Example of DE with TMM normalization (See CompCodeR Vignette)
runDiffExp(data.file = "B_625_625_5spc_repl1.rds", 
           result.extent = "voom.limma", Rmdfunction = "voom.limma.createRmd", 
           output.directory = ".", norm.method = "TMM")
runDiffExp(data.file = "B_625_625_5spc_repl1.rds", 
           result.extent = "edgeR.exact", Rmdfunction = "edgeR.exact.createRmd", 
           output.directory = ".", norm.method = "TMM", 
           trend.method = "movingave", disp.type = "tagwise")
runDiffExp(data.file = "B_625_625_5spc_repl1.rds", result.extent = "ttest", 
           Rmdfunction = "ttest.createRmd", 
           output.directory = ".", norm.method = "TMM")

B_625_625_5spc_repl1 = readRDS('B_625_625_5spc_repl1.rds')
count = data.frame(B_625_625_5spc_repl1@count.matrix)

data = as_tibble(t(count),rownames=NA)
data = data%>%
  mutate(ID=row_number()) %>%
  mutate(condition = c(rep('condition_1',5),rep('condition_2',5)))
  
metadata = data %>%
  select(ID,condition)

set.seed(1246)

#https://academic.oup.com/bioinformatics/article/29/22/2877/313226?login=true

#****************************************************************************
# Plot ------------------------------------------------------------------------
umap_fit = data %>%
  select(where(is.numeric))%>%
  column_to_rownames('ID')%>%
  scale()%>%
  umap(n_neighbors=5)

umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(metadata, by="ID")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = condition,))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")
ggsave("UMAP_Normal.png")

df = data %>%
  select(-ID,-condition)
pca_res <- prcomp(df, scale. = TRUE)
autoplot(pca_res, data = data, colour = 'condition')
ggsave("PCA_Normal.png")

#****************************************************************************  
# quantile normalization -------------------------------------------------------
#Developed for gene expression microarrays but applied to many high-throughput
#techniques including rna-seq and proteomics
#https://www.nature.com/articles/s41598-020-72664-6

#The quantile normalization (QN) procedure is simple (Fig. 1A): it involves 
#first ranking the gene of each sample by magnitude, calculating the average 
#value for genes occupying the same rank, and then substituting the values of 
#all genes occupying that particular rank with this average value. The next step
#is to reorder the genes of each sample in their original order. 


#All- Naive application of QN


#perform quantile normalization
df_norm <- as.data.frame(normalize.quantiles(as.matrix(df)))

pca_res <- prcomp(df_norm, scale. = TRUE)
autoplot(pca_res, data = data, colour = 'condition')
ggsave("PCA_All_Quantile.png")

umap_fit = df_norm %>%
  scale()%>%
  umap(n_neighbors=5)

umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(metadata, by="ID")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = condition,))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")
ggsave("UMAP_All_Quantile.png")

#Class-specific - QN applied on each subclass
df_1 = data %>%
  filter(condition == 'condition_1') %>%
  select(-ID,-condition)

df_norm_1 <- as.data.frame(normalize.quantiles(as.matrix(df_1)))

df_2 = data%>%
  filter(condition=='condition_2')%>%
  select(-ID,-condition)

df_norm_2 <- as.data.frame(normalize.quantiles(as.matrix(df_2)))

df_norm = rbind(df_norm_1,df_norm_2)

pca_res <- prcomp(df_norm, scale. = TRUE)
autoplot(pca_res, data = data, colour = 'condition')
ggsave("PCA_Class_Quantile.png")

umap_fit = df_norm %>%
  scale()%>%
  umap(n_neighbors=5)

umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(metadata, by="ID")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = condition,))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")
ggsave("UMAP_Class_Quantile.png")


#****************************************************************************
# upper quantile---------------------------------------------------------------

normalization_factor <- as_tibble(calcNormFactors(count,method="upperquartile"))
df_norm = mapply(`*`,df,normalization_factor)

pca_res <- prcomp(df_norm, scale. = TRUE)
autoplot(pca_res, data = data, colour = 'condition')
ggsave("PCA_Upper_Quantile.png")

umap_fit = df_norm %>%
  scale()%>%
  umap(n_neighbors=5)

umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(metadata, by="ID")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = condition,))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")
ggsave("UMAP_Upper_Quantile.png")


#****************************************************************************
# z-score----------------------------------------------------------------------


#****************************************************************************
# min-max----------------------------------------------------------------------



#****************************************************************************
# log2-------------------------------------------------------------------------



#****************************************************************************
#trimmed or whitening normalization

#****************************************************************************
# Other -----------------------------------------------------------------------
#Combat and other normalization technique can be considered normalization tech-
#niques
#GC content normalization
#Omics specific (CPM, TPMM,RPKM) for rna seq

#****************************************************************************
#Evaluation using gPCA t-test (alpha=0.05)
#Extension of the approach using different techniques
#edgeR calling and comparison




#guide deseq edgeR
#https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25

#****************************************************************************¸
# RNA-Seq Normalization--------------------------------------------------------
#within sample RPKM - QUANTILE
#Between sample CMP,TPMM, upper quantile -- great caution if ML is involved

#Visualization
#PCA (principal component analysis)
#MDS (multidimensional scaling)#
#Boxplot, Density Plots, Heatmaps

#*******************************************************************************
# RNA Seq -----------------------------------------------------------------
#*******************************************************************************
## Import data in RStudio
counts <- read.table("Data/pbmc_count.txt", sep = "\t", header = TRUE)
counts1 <- counts[,-1]
rownames(counts1) <- sapply(strsplit(as.character(counts$Ensembl),"[.]"),function(x)(x[1]))




## Annotation
# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# ensembl <- as.character(rownames(dt))
# ann <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), 
#              filters = c("ensembl_gene_id"),values = ensembl, mart = mart)
# write.table(ann, file = "annotation_hsa.txt", sep = "\t", row.names = F, col.names = T)

ann <- read.table("Data/annotation_hsa.txt", sep = "\t", header = TRUE)
ann <- merge(ann,counts1,by.x = 1,by.y = 0)
ann <- ann[!duplicated(ann$ensembl_gene_id),]

# Barplot Library size
dt <- ann[,-c(1:2)]
rownames(dt) <- ann$ensembl_gene_id

lz <- data.frame(apply(dt,2,sum))
colnames(lz) <- "library_size"
bar <- ggplot(lz, aes(x=rownames(lz), y=library_size)) +
  geom_bar(stat="identity")+theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="")
print(bar)

# Metadata
NP <- c(rep("N",3),rep("P",3))
groups <- factor(NP, levels = c("N","P"), labels = c("N","P"))

## Filtering the data (filter lowly expressed genes)
th <- min(table(groups)) - 1
data <- as.matrix(dt)
rec <- as.matrix(data[rowSums(cpm(data) > 1) > th, ]) 
data <- rec
print(dim(data))

## TMM Normalization
normdata <- tmm(data) 
ann <- merge(ann[,c(1:2)], as.data.frame(normdata),by.x = 1, by.y = 0)

df <-  melt(normdata, variable_name = 1)[,2:3]
colnames(df) <- c("sample", "TMM")
df$condition <- c(rep("N", length(grep("N",df$sample))),rep("P", length(grep("P",df$sample))))
box <- ggboxplot(df, "sample", "log2(TMM)", color = "condition",palette = c("#00AFBB", "#E7B800"))
box

## Visualization of the samples: PCA
pca_cl <- prcomp(t(log2(normdata+1)))
df_out <- as.data.frame(pca_cl$x)

name <- row.names(df_out)
shape <- NP

percentage <- round(pca_cl$sdev/sum(pca_cl$sdev) * 100, 2)
percentage <- paste(colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep=""))

p <- ggplot(df_out,aes(x=PC1,y=PC2,colour=groups,shape = shape)) +
  geom_text_repel(aes(label=name), size = 3)+ 
  xlab(percentage[1]) + ylab(percentage[2])
p <- p+geom_point(aes(colour=groups), size=3, stroke = 1) + ggtitle("PCA Plot") + 
  labs(colour = "",shape = "") +
  theme_bw()
print(p)

## DE
# Make design matrix
design <- model.matrix(~0 + groups)
rownames(design) <- colnames(data)
colnames(design) <- gsub("groups", "", colnames(design))

# For easy manipulation, we put the data into a DGEList object.
y <- DGEList(counts=data)

# TMM normalization is applied to this dataset to account for compositional difference between the libraries.
y <- calcNormFactors(y)
y$samples

# Estimating the dispersion
y <- estimateGLMRobustDisp(y,design)

# Fit the model
fit <- glmFit(y,design)
lrt <- glmLRT(fit, contrast = c(-1,1))

# Total number of differentially expressed genes at 5% FDR
print(summary(decideTests(lrt)))

# Plot log-fold change against log-counts per million, with DE genes highlighted.
plotMD(lrt)
abline(h=c(-1, 1), col="blue") # The blue lines indicate 2-fold changes.
nr <- dim(data)[1]

res <- topTags(lrt, n = nr)$table
ann <- merge(ann, res, by.x = 1, by.y = 0)
ann <- ann[order(ann$FDR),]
sum(ann$FDR < 0.05)

# Save data
write.table(ann, file = "ResultsDE.txt", sep = "\t", row.names = F, col.names = T)