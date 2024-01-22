##WGCNA
install.packages('WGCNA')
install.packages("gridExtra")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")
BiocManager::install("preprocessCore")
BiocManager::install("GEOquery")
devtools::install_github("kevinblighe/CorLevelPlot")

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)


#get annotation - i used the 110 bulid (Release 44 (GRCh38.p14))
hub = AnnotationHub()                  
ensdb_query <- query(hub, c("EnsDb", "sapiens", "110"))
ensdb_query

#retreieve record
ensdb_110 = ensdb_query[["AH113665"]]

# Extract transcript and gene information
tx_data <- transcripts(ensdb_110, return.type = "DataFrame")


# Create the tx2gene data.frame
tx2gene <- tx_data[, c("tx_id", "gene_id")]

#View
view(data.frame(tx2gene))

#prepare data
quants_dir = "Data/"

#get  quant file path
quant_files = list.files(quants_dir, pattern = "quant.sf$", recursive = TRUE, full.names = TRUE)

#view path
quant_files

#get file directory path
quant_dirs <- list.files(quants_dir, pattern = ".salmon$", full.names = TRUE)

#get sample names - match with what is in metadata
#sample_names <- gsub("(.*)_S[0-9]+\\.salmon$", "\\1", basename(quant_dirs))

#OR

sample_names <- gsub("_.*\\.salmon$", "", basename(quant_dirs))

#combine sample names to the quant file path into a named list
names(quant_files) <- sample_names
quant_files

#pass to txi
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene,ignoreTxVersion = TRUE)

#import metdadata csv
# Specify the path to your CSV file
file_path <- "Meta/Final_meta_tissue_healvsamp.csv"

# Read the CSV file into a data frame
meta <- read.csv(file_path, row.names = 1)

View(meta)

##subset the metadata  to only those variables of interest
head(meta)
meta_sub <- meta[,c(1:5,8:12,46:47,111)] #select desired columns
head(meta_sub)

#check if colnmaes match
all(colnames(txi$counts) %in% rownames(meta_sub))
all(colnames(txi$counts) == rownames(meta_sub))

#create DEseq2 object
dds <- DESeqDataSetFromTximport(txi, colData = meta_sub, design = ~1)
dds
View(counts(dds))


##detect outlier
data = (counts(dds))
#any outlier gene or samples?
gsg <- goodSamplesGenes(t(data))

summary(gsg)
gsg$allOK #if true, then all genes and samples passed. if false, then either gene or sample( or both) failed

#shows how many genes/samples are true or false for outlier
table(gsg$goodGenes)
table(gsg$goodSamples)

#remove genes detected as outliers
dds <- dds[gsg$goodGenes == TRUE,]
dds

# detect outlier samples - hierarchical clustering - method 1
##detect outlier
data = (counts(dds))
nrow(data)

htree <- hclust(dist(t(data)), method = "average")
plot(htree) ##you can remove samples not well clustered

# pca - method 2

#calculate PC
pca <- prcomp(t(data))
pca.dat <- pca$x

#calculate vriance explained by each PC
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

#convert the PC compoenets to dataframe
pca.dat <- as.data.frame(pca.dat)

#use ggplot2 to plot
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# exclude outlier samples from the data
samples.to.be.excluded <- c('LKRNA003-76', 'LKRNA003-82', 'LKRNA003-57', 'LKRNA003-86')
dds_sub <- dds[,!(colnames(dds) %in% samples.to.be.excluded)]

#counts2 = (counts(dds_sub))

# exclude outlier samples from the metadata
colData <- meta_sub %>% 
  dplyr::filter(!row.names(.) %in% samples.to.be.excluded)  

### further filtering - remove all genes with counts < 15 in more than 75% of samples (42*0.75=31.5)
## suggested by WGCNA on RNAseq FAQ
dds75 <- dds_sub[rowSums(counts(dds_sub) >= 15) >= 32,]
view(counts(dds75)) 
nrow(dds75) 

#OR

##prefiltering - remove rows where total sum is less than 10
#dds10_2 <- dds_subset[rowSums(DESeq2::counts(dds_subset)) > 10]
#dds10_2
#nrow(dds10)

# perform variance stabilization
dds_norm <- vst(dds75)
view(counts(dds_norm))

#vew plot
plotPCA(dds_norm, intgroup="Batch", pcsToUse = 1:2)

#remove batch effect
mat <- assay(dds_norm)
mm <- model.matrix(~1, colData(dds_norm))
mat <- limma::removeBatchEffect(mat, batch=dds_norm$Batch, design=mm)
assay(dds_norm) <- mat

#vew plot
plotPCA(dds_norm, intgroup="Batch", pcsToUse = 3:4)

# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

###WGCNA - NETWORK construction
# Choose a set of soft-thresholding powers
#This uses default of power  = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
sft <- pickSoftThreshold(norm.counts,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)

?pickSoftThreshold
##or CHOOSE from predefined power cutoff
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) + 
  geom_point() + # Plot the points
  geom_text(nudge_y = 0.1) + # We'll put the Power labels slightly above the data points
  geom_hline(yintercept = 0.8, color = 'red') + # We will plot what WGCNA recommends as an R^2 cutoff
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

a1/a2

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 9 #picked power
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 11000, #specify how many genes in a block- here she used enoughh size to contain all genes in one block 
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,#color instead of numbers
                          randomSeed = 1234,
                          verbose = 3)

cor <- temp_cor 

?blockwiseModules
# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# grey module = all genes that doesn't fall into any modules were assigned to the grey module
#the unmerged has more colors than the merged which indicates that some similar modules were merged. We will use the merged module for further downstream analysis

#Biological questions
#what are the genes or clusters of genes (modules) significanty asociated withh COVID-19 individuals?
#Identify genes that are significantly associated with severe COVID-19 cases


#Relate/associate modules to trait - which module have significant association with trait of interest (outcome for me, COVID-19 and severity for above)

#binarize the categorical variable of interest
traits <- colData %>% 
  mutate(outcome_bin = ifelse(grepl('Healed', outcome), 1, 0)) %>%  #if amputated (value of interest in this case), give a value of 1, or else (healed) give a value of 0
  dplyr::select(11,14) #select column of interest ie outcome_bin


#OR - better if more than one category, here columns are made for each variable and column with a varaible for that column is assigned 1 while the rest is assigned 0
#example - column healthy will have value of 1 if the sample is healthy and O value for other samples not healty (either seevre, ICU etc)

# binarize categorical variables

#colData$severity <- factor(colData$severity, levels = c("Healthy", "Convalescent", "ICU", "Moderate", "Severe"))

#severity.out <- binarizeCategoricalColumns(colData$severity,
#                                        includePairwise = FALSE,
#                                          includeLevelVsAll = TRUE,
#                                         minCount = 1)


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

#calculate corrrelation between module eigenegene and trait(s) of interest using pearson correlaton
module.trait.corr <- cor(module_eigengenes, traits, use = 'p')

#calaculate pvalues for these correlations
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names') #convert the column names to rownames if not converted


CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[18:22], #trait columns
             y = names(heatmap.data)[1:17], #Module columns
             col = c("blue1", "skyblue", "white", "pink", "red")) #based on the number of trait columns

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[17:18], #trait columns
             y = names(heatmap.data)[1:16], #Module columns
             col = c("blue1", "skyblue")) 


#level of significacne is determined by the astericks, 3 astericks shows high significance to that trait of interest (COVID-19 patients or Amputated)
#extract genes in modules with high significane (3 astericks)

#extract the genes that are part of the significant module , in this case turqoise 
module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()


#The result above will give you same thing as the limma fxn analysis below

# 6B. Intramodular analysis: Identifying driver genes ---------------
# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.


#calculate corrrelation between module eigenegene and genes expressio profile using pearson correlaton
module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')


#calaculate pvalues for this correlation/measures
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

# Using the module membership measures you can identify genes with significantly high module membership in interesting modules by looking at the p-values
#you can extract these genes and look further into it


#genes associated with trait of interest ie genes that have a high significance for trait of interest 

#correlate expression data with trait of interest
gene.signf.corr <- cor(norm.counts, traits$outcome_bin, use = 'p')

#calculate pvalue of correlation
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

#top 25 genes significantly associated with Healed patients
top_25_healed <- gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  filter(V1 < 0.05) %>%
  head(25)

nrow(top_25_healed)
#box

subset_sig_box = t(norm.counts)[c("ENSG00000189143"), ]
subset_sig_box = data.frame(subset_sig_box)%>%
  rownames_to_column(var = "sample_name")%>%
  rename(Expression = subset_sig_box) %>% #rename subset_sig_box to Expression
  mutate(SampleID = str_replace_all(sample_name, "-", ".")) %>% 
  dplyr::left_join (dplyr::select(colData_df,sample_name,Type,outcome), by = "sample_name") #only add type and outcome to the merged dataframe



ggplot(subset_sig_box, aes(x =outcome, y = Expression, fill = outcome)) +
  geom_boxplot() +
  labs(title = "Gene - DPM1",
       x = "Grouping Variable",
       y = "Expression level") +
  theme_minimal()

##limma
module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)


## Which modules have biggest differences across treatment groups?

#We can also see if our eigengenes relate to our metadata labels. 
#First we double check that our samples are still in order.
all.equal(rownames(colData), rownames(module_eigengenes))

# Create the design matrix from the `outcome` variable
des_mat <- model.matrix(~ colData$outcome)


#Run linear model on each module.
#Limma wants our tests to be per row, so we also need to transpose so the eigengenes are rows


# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)


#Apply multiple testing correction and obtain stats in a data frame. 

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")


#Let's take a look at the results. 
#They are sorted with the most significant results at the top.

head(stats_df)


#Module midnightblue seems to be the most differentially expressed across `time_point` groups. 
#Now we can do some investigation into this module. 

## Let's make plot of module midnightblue

#As a sanity check, let's use `ggplot` to see what module midnightblue eigengene looks like between treatment groups. 

#First we need to set up the module eigengene for this module with the sample metadata labels we need. 

colData_df = colData %>%
rownames_to_column("sample_name")

module_midnight_df <- module_eigengenes %>%
  rownames_to_column("sample_name") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(colData_df %>%
                      dplyr::select(sample_name, outcome, wagner_grade),
                    by = "sample_name")

#Now we are ready for plotting. 

ggplot(module_midnight_df,aes(x = outcome,y = MEmidnightblue,color = outcome)) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()

## What genes are a part of module midnight?

#If you want to know which of your genes make up a modules, you can look at the `$colors` slot. 
#This is a named list which associates the genes with the module they are a part of. 
#We can turn this into a data frame for handy use. 

gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))


#Now we can find what genes are a part of module 19. 
gene_module_midnight <- gene_module_key %>%
  dplyr::filter(module == "MEmidnightblue")

## Make a custom heatmap function

#We will make a heatmap that summarizes our differentially expressed module.
#Because we will make a couple of these, it makes sense to make a custom function for making this heatmap. 

make_module_heatmap <- function(module_name,
                                expression_mat = norm.counts,
                                metadata_df = colData_df,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts
  # metadata_df: a data frame with refinebio_accession_code and time_point
  #              as columns. Default is `metadata'
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key'
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes'
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its refinebio_accession_code
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("sample_name")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(outcome, sample_name, wagner_grade, Patient) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "sample_name") %>%
    # Arrange by patient and time point
    dplyr::arrange(outcome, Patient) %>%
    # Store sample
    tibble::column_to_rownames("sample_name")
  
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    outcome = col_annot_df$outcome,
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in time_point
    col = list(outcome = c("Amputated" = "#f1a340", "Healed" = "#998ec3"))
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = FALSE,
                                     show_column_names = FALSE
  )
  
  # Return heatmap
  return(heatmap)
}

## Make module heatmaps
#Let's try out the custom heatmap function with module 19 (our most differentially expressed module).
#we can also plot heatmap of other non differentially expressed modules by changing the module name

mod_midnight_heatmap <- make_module_heatmap(module_name = "MEblue")

# Print out the plot
mod_midnight_heatmap

#From the barplot portion of our plot, we can see healed samples tend to have higher expression values for the module midnight eigengene. 
#In the heatmap portion, we can see how the individual genes that make up module midnight are overall higher in healed than in the Amputated samples

