#Description and effect of normalization techniques
library(tidyverse)
library(ggfortify)
library(preprocessCore)
library(edgeR)
library(caret)

#DATASET CompCODEr
browseVignettes("compcodeR")

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
# variance equal to one ^40 .

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
                                   dispersions ='auto', 
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

df = data %>%
  select(-ID,-condition)

set.seed(1246)

#https://academic.oup.com/bioinformatics/article/29/22/2877/313226?login=true

#****************************************************************************
# Plot ------------------------------------------------------------------------
umap_fit = data %>%
  select(where(is.numeric))%>%
  column_to_rownames('ID')%>%
  #scale()%>%
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

pca_res <- prcomp(df, scale. = F)
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
df_t <- t(df)

#All- Naive application of QN
#perform quantile normalization
df_norm_t <- as.data.frame(normalize.quantiles(as.matrix(df_t)))
df_norm <- t(df_norm_t)

pca_res <- prcomp(df_norm, scale. = F)
autoplot(pca_res, data = data, colour = 'condition')
ggsave("PCA_All_Quantile.png")

umap_fit = df_norm %>%
  #scale()%>%
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

df_1_t <- t(df_1)
df_norm_1_t <- as.data.frame(normalize.quantiles(as.matrix(df_1_t)))
df_norm_1 <- t(df_norm_1_t)

df_2 = data%>%
  filter(condition=='condition_2')%>%
  select(-ID,-condition)

df_2_t = t(df_2)
df_norm_2_t <- as.data.frame(normalize.quantiles(as.matrix(df_2_t)))
df_norm_2 <- t(df_norm_2_t)

df_norm = rbind(df_norm_1,df_norm_2)

pca_res <- prcomp(df_norm, scale. = F)
autoplot(pca_res, data = data, colour = 'condition')
ggsave("PCA_Class_Quantile.png")

umap_fit = df_norm %>%
  #scale()%>%
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

pca_res <- prcomp(df_norm, scale. = F)
autoplot(pca_res, data = data, colour = 'condition')
ggsave("PCA_Upper_Quantile.png")

umap_fit = df_norm %>%
  #scale()%>%
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

#Per class upper quantile
df_1 =  data %>%
  filter(condition == 'condition_1') %>%
  select(-ID,-condition)


normalization_factor_1 <- as_tibble(calcNormFactors(t(df_1),method="upperquartile"))
df_norm_1 = mapply(`*`,df_1,normalization_factor_1)


df_2 = data%>%
  filter(condition=='condition_2')%>%
  select(-ID,-condition)

normalization_factor_2 <- as_tibble(calcNormFactors(t(df_2),method="upperquartile"))
df_norm_2 = mapply(`*`,df_2,normalization_factor_2)

df_norm = rbind(df_norm_1,df_norm_2)

pca_res <- prcomp(df_norm, scale. = F)
autoplot(pca_res, data = data, colour = 'condition')
ggsave("PCA_Upper_Quantile_Class.png")


umap_fit = df_norm %>%
  #scale()%>%
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
ggsave("UMAP_Upper_Quantile_Class.png")

#****************************************************************************
# z-score----------------------------------------------------------------------
#Normalize each feature computing z = (x-MEAN(x))/SD(x)
df_norm = sapply(df, function(df) (df-mean(df))/sd(df))
pca_res <- prcomp(df_norm, scale. = F)
autoplot(pca_res, data = data, colour = 'condition')
ggsave("PCA_Z_Score.png")


umap_fit = df_norm %>%
  #scale()%>%
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
ggsave("UMAP_Z_Score.png")

#****************************************************************************
# min-max----------------------------------------------------------------------
process = preProcess(df, method=c("range"))
df_norm <- predict(process, df)

pca_res <- prcomp(df_norm, scale. = F)
autoplot(pca_res, data = data, colour = 'condition')
ggsave("PCA_Min_Max.png")


umap_fit = df_norm %>%
  #scale()%>%
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
ggsave("UMAP_Min_Max.png")


#****************************************************************************
# log2-------------------------------------------------------------------------
df_norm = sapply(df, function(df)(log(df+1)))
pca_res <- prcomp(df_norm, scale. = F)
autoplot(pca_res, data = data, colour = 'condition')
ggsave("PCA_LOG.png")

umap_fit = df_norm %>%
  #scale()%>%
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
ggsave("UMAP_Log.png")
