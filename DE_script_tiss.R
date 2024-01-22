if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("ensembldb")
BiocManager::install("AnnotationHub")
BiocManager::install("DESeq2")

install.packages("tidyverse")
install.packages("pheatmap")
install.packages("patchwork")


install.packages("ggrepel")
library(tximport)
library(ensembldb)
library(AnnotationHub)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(patchwork)
library(ggrepel)


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
data.frame(tx2gene)

#prepare data
quants_dir = "Data/"

#get  quant file path
quant_files = list.files(quants_dir, pattern = "quant.sf$", recursive = TRUE, full.names = TRUE)

#view path
quant_files

#get file directory path
quant_dirs <- list.files(quants_dir, pattern = ".salmon$", full.names = TRUE)

#get sample names - match with what is in metadata
sample_names <- gsub("(.*)_S[0-9]+\\.salmon$", "\\1", basename(quant_dirs))

#OR

sample_names <- gsub("_.*\\.salmon$", "", basename(quant_dirs))

#combine sample names to the quant file path into a named list
names(quant_files) <- sample_names
quant_files

#pass to tximport
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene,ignoreTxVersion = TRUE)

colnames(txi$counts)

txi$counts



#OR 

#txi_2 <- tximport(quant_files, type="salmon", tx2gene=tx2gene, countsFromAbundance="lengthScaledTPM", ignoreTxVersion = TRUE)

#view 
txi.dat = data.frame(txi)

#write.table(txi.dat, file="Data/txi.txt", sep="\t", quote=F, col.names=NA)

#import metdadata csv
# Specify the path to your CSV file
file_path <- "Meta/Final_meta_tissue_healvsamp.csv"

# Read the CSV file into a data frame
meta <- read.csv(file_path, row.names = 1)

View(meta)
#################################
#make sure column names match
all(colnames(txi$counts) %in% rownames(meta))
all(colnames(txi$counts) == rownames(meta))

#reorder to match
idx = match(colnames(txi$counts), rownames(meta))
re_meta = meta[idx, ]
View(re_meta)

all(colnames(txi$counts) == rownames(re_meta))

##Prefiltering
#########check DESEQ data in old laptop

#create deseq2 object
dds_norm <- DESeqDataSetFromTximport(txi, colData = meta, design = ~Batch + outcome)
dds_norm

View(counts(dds))



##prefiltering - remove rows (genes) with only 0 or 1 read
#dds <- dds[ rowSums(counts(dds)) > 1]
#dds

#or filter rows(genes) where total sum is less than 10
dds_norm <- dds_norm[rowSums(DESeq2::counts(dds_norm)) > 10]
dds_norm

View(counts(dds_norm))

##perform normalization
dds_norm <- estimateSizeFactors(dds_norm)
#sizeFactors(dds)
normalizationFactors(dds_norm)

#retrieve normalized count
#normalized_counts <- counts(dds, normalized=TRUE)
#View(normalized_counts)
#write.table(normalized_counts, file="Data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)
##Deseq2 doesn't use normalized count, however they can be useful for downstram visualization


#Sample level QC - PCA and heatmap
##These unsupervised clustering methods are run using log2 transformed normalized counts (regularized) or VST
#The DESeq2 vignette suggests large datasets (100s of samples) to use the variance-stabilizing transformation (vst) instead of rlog for transformation of the counts, 
#since the rlog function might take too long to run and the vst() function is faster with similar properties to rlog.

#rld <- rlog(dds, blind=TRUE) #if less samples

#OR VST - blind = true since it is quality assesment

vsdata <- vst(dds_norm, blind = TRUE)

#remove batch effect
mat <- assay(vsdata)
mm <- model.matrix(~outcome, colData(vsdata))
mat <- limma::removeBatchEffect(mat, batch=vsdata$Batch, design=mm)
assay(vsdata) <- mat


#Plot PCA
plotPCA(vsdata, intgroup="Batch", pcsToUse = 1:2)
plotPCA(vsdata, intgroup="Batch", pcsToUse = 3:4)

A = plotPCA(vsdata, intgroup="outcome", pcsToUse = 1:2)
B = plotPCA(vsdata, intgroup="outcome", pcsToUse = 3:4)

C = A + ggtitle('P1/P2') + labs(color = "outcome") + stat_ellipse()
D = B + ggtitle('P3/P4') + labs(color = "outcome") + stat_ellipse()

C + D + plot_annotation(title = "PCA plot for Outcome - Tissue") 

#A + ggtitle('P1/P2') + labs(color = "Batch") + stat_ellipse()

#heatmap
library("pheatmap")

#select the top 20 genes based on the mean expression
select <- order(rowMeans(counts(dds_norm,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_norm)[,c("Type","outcome", "Batch")])
pheatmap(assay(vsdata)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


#sample to sample distance
library("RColorBrewer")
sampleDists <- dist(t(assay(vsdata)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsdata$outcome, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



###Add labels to plot
library(ggplot2)

# Assuming "Sample" is the column containing sample names in vsdata - change intgroup and color to Type if you want to see sample type relationship
sample_names <- vsdata$Patient

# Generate PCA plot - also check pc2 and 3 for outlier identification
pca_plot <- plotPCA(vsdata, intgroup = "outcome", pcsToUse = 1:2)

# Extract data from the PCA plot
pca_data <- as.data.frame(pca_plot$data)

# Combine PCA data with sample names
pca_data_with_names <- cbind(pca_data, Sample = sample_names)

# Create a custom PCA plot with sample names
custom_pca_plot <- ggplot(pca_data_with_names, aes(x = PC1, y = PC2, color = outcome)) +
  geom_point() +
  geom_text(aes(label = Sample), vjust = 1, hjust = 1, size=3) +  # Add sample names
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  stat_ellipse() +
  theme_minimal()

# Display the custom PCA plot
print(custom_pca_plot)


#custimize PCA plot
pcaData <- plotPCA(vsdata, intgroup=c("outcome", "Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Batch, shape=outcome)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

###HEATMAP ... will cluster like the PCA above

#compute pairwise correlation values
vst_mat <- assay(vsdata) #create matrix
vst_cor = cor(vst_mat)

#lok at columns and row names, #do the names match?
head(vst_cor)
head(meta)

#plot pheatmap
annotation_col = data.frame(
  `Outcome` = as.factor(meta$outcome),
  `Type` = as.factor(meta$Type),
  `Wagner` = as.factor(meta$wagner_grade),
  `Batch` = as.factor(meta$Batch),
  check.names = FALSE
)
rownames(annotation_col) = rownames(meta)

heat.colors <- RColorBrewer::brewer.pal(6, "YlOrRd")
pheatmap(vst_cor, annotation = annotation_col, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20)



###DEANAlsyis
## Create DESeq2Dataset object - fix design based on results above
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ Batch + outcome)

#prefilter
dds <- dds[rowSums(DESeq2::counts(dds)) > 10]
dds

#run Deseq2
dds = DESeq(dds) #run deseq2

#total number of raw counts persample
colSums(counts(dds))

#totalnumber of normalized counts per sample
colSums(counts(dds, normalized=T))

#create dispersion plot to check if model fits right
plotDispEsts(dds)

#plot PCA of result
vsdata = vst(dds, blind=FALSE)


A = plotPCA(vsdata, intgroup="outcome", pcsToUse = 1:2)
B = plotPCA(vsdata, intgroup="Batch", pcsToUse = 1:2)

A | B


#show results = baseline is healed - lower or higher in healed compared to amputated
#Positive log2 fold changes indicate genes upregulated/down in "Healed" compared to "amputated," 
res = results(dds, contrast = c('outcome', 'Healed', 'Amputated'))

##Apply LFC shrink
resLFC <- lfcShrink(dds, coef="outcome_Healed_vs_Amputated", type="apeglm", res = res)

?lfcShrink
resultsNames(dds) #run this to know the correct coef value to use, it is st

#plot MA plot using shrunken result
plotMA(resLFC, ylim=c(-2,2))

#check both
summary(resLFC)
summary(res)


#explore result column
mcols(res)
head(res, n=10)
summary(res)




##Since I have fewer genes, I will only select significant ones based on the p-value
# Create a tibble of signifciant genes results
#sigs2_tb <- sigs_2 %>%
#  data.frame() %>%
 # rownames_to_column(var="gene") %>% 
 # as_tibble()


library(annotables)
library(dplyr)
grch38

#select only few desired column of grch38
grch = grch38[, c("ensgene", "symbol", "description")] 
grch_f = distinct(grch)
  
nrow(grch38)
nrow(grch_f)

all_res <- data.frame(resLFC) %>%
  rownames_to_column(var="ensgene") %>% 
  left_join(
    y = grch_f[, c("ensgene", "symbol", "description")],
    by = "ensgene")
nrow(all_res)

##extract filtered result  by padj and/or basemean then arrange by padj
all_res_sig = subset(all_res, padj < 0.05)

all_res_sig = all_res_sig %>%
  arrange(padj)

View(all_res_sig)

#subset normalzed counts to significant genes
normalized_counts <- counts(dds, normalized=T)

#plot the heatmap to be plotted by the gene name instead
sig_norm_counts = normalized_counts[all_res_sig$ensgene, ]
sig_norm_counts_2= sig_norm_counts[1:50, ]

label = all_res_sig$symbol[1:50]

##expression heatmap
#library(RColorBrewer)
heat_colors <- brewer.pal(6, "YlOrRd")
#display.brewer.all()
?pheatmap
#run heatmap
pheatmap(sig_norm_counts_2, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = T,
         annotation = select(meta, c("outcome", "Type", "wagner_grade", "Batch")),
         #border_color = NA, 
         #fontsize = 10, 
         scale = "row", 
         #fontsize_row = 10, 
         labels_row = label,
         height = 30)

#Volcano Plot
all_res_sig_2 <- all_res_sig %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)

#all_res_sig_3 <- all_res_sig %>% 
  #mutate(threshold = padj < 0.05 & log2FoldChange >= 1.5)

summary(res)
?results
#volcano plot
ggplot(all_res_sig_2) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle("HealedvsAmputated - Tissue") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
   theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  


nrow(all_res_sig)
nrow(all_res_sig_2)

#label the top 20 significant genes on the volcano plot
## Create an empty column to indicate which genes to label
all_res_tb  <- all_res_sig  %>% dplyr::mutate(genelabels = "")
all_res_tb$genelabels[1:15] <- as.character(all_res_tb$symbol[1:15])
all_res_tb<- all_res_tb %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)

ggplot(all_res_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("HealedvsAmputated - Tissue") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

#view signicant results
View(all_res_sig)

#view normalized counts for signifcant genes
View(sig_norm_counts)
ncol(sig_norm_counts)


all_res_tb  <- all_res_sig  %>% dplyr::mutate(genelabels = "")
all_res_tb$genelabels[1:15] <- as.character(all_res_tb$symbol[1:15])

#Visulaize top 20 genes - make sure 2:46 columns includes the last column
top_20 = data.frame(sig_norm_counts)[1:20, ] %>%
  rownames_to_column(var = "gene")
top_20 = gather(top_20, key = "SampleID", value = "normalized_counts", 2:47)
top_20


#convert meta rowname to a column named sampleID and then innerjoin
meta_wt = rownames_to_column(meta, var = "SampleID") 
meta_wt$SampleID = gsub("-", ".", meta_wt$SampleID)

top_20 = inner_join(top_20, meta_wt, by = "SampleID")

#plot - you can use x = gene if you want the gene symbol
ggplot(top_20) +
  geom_point(aes(x = gene, y = normalized_counts, color = outcome)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))


write.table(all_res_sig, file="Data/all_res_sig.txt", sep="\t", quote=F, col.names=NA)


###Using DEGreport to summarize results _ recheck options

DEGreport::degPlot(dds = dds, res = resLFC, n = 20, xs = "outcome", group = "outcome") # dds object is output from DESeq2

DEGreport::degVolcano(
  data.frame(res[,c("log2FoldChange","padj")]), # table - 2 columns
  plot_text = data.frame(res[1:10,c("log2FoldChange","padj")])) # table to add names

# Available in the newer version for R 3.4
DEGreport::degPlotWide(counts = dds, genes = row.names(res)[1:20], group = "outcome")



#visualize only the genes listed here

top_20 = data.frame(sig_norm_counts)[1:20, ] %>%
  rownames_to_column(var = "gene")
top_20 = gather(top_20, key = "SampleID", value = "normalized_counts", 2:47)
top_20

g_interest = data.frame(sig_norm_counts) %>%
            rownames_to_column(var = "gene")

row_interest <- g_interest[g_interest$gene %in% c("ENSG00000096696", "ENSG00000158055", "ENSG00000143631", "ENSG00000133710", "ENSG00000203782",
                                                  "ENSG00000205420", "ENSG00000118898", "ENSG00000189143", "ENSG00000105141", "ENSG00000128422",
                                                  "ENSG00000186832", "ENSG00000124882", "ENSG00000167880", "ENSG00000167768", "ENSG00000164687",
                                                  "ENSG00000163207", "ENSG00000196586", "ENSG00000186847", "ENSG00000179148"), ]


 row_interest = gather(row_interest, key ="SampleID", value = "normalized_counts", 2:47)
 row_interest

 #convert meta rowname to a column named sampleID and then innerjoin
 meta_wt = rownames_to_column(meta, var = "SampleID") 
 meta_wt$SampleID = gsub("-", ".", meta_wt$SampleID)
 
 row_interest = inner_join(row_interest, meta_wt, by = "SampleID")
 
 #plot - you can use x = gene if you want the gene symbol
 ggplot(row_interest) +
   geom_point(aes(x = gene, y = normalized_counts, color = outcome)) +
   scale_y_log10() +
   xlab("Genes") +
   ylab("log10 Normalized Counts") +
   ggtitle("Epidermis development Genes") +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
   theme(plot.title = element_text(hjust = 0.5))


 #plot the heatmap to be plotted by the gene name instead
 sig_norm_counts = normalized_counts[all_res_sig$ensgene, ]
 sig_norm_counts_2= sig_norm_counts[1:50, ]
 
 subset_sig_norm = sig_norm_counts[c("ENSG00000096696", "ENSG00000158055", "ENSG00000143631", "ENSG00000133710", "ENSG00000203782",
                                     "ENSG00000205420", "ENSG00000118898", "ENSG00000189143", "ENSG00000105141", "ENSG00000128422",
                                     "ENSG00000186832", "ENSG00000124882", "ENSG00000167880", "ENSG00000167768", "ENSG00000164687",
                                     "ENSG00000163207", "ENSG00000196586", "ENSG00000186847", "ENSG00000179148"), ]

 
 label = all_res_sig$symbol[1:50]
 
 class(sig_norm_counts)
 
 any(is.na(heat_interest))
 any(is.infinite(heat_interest))
 
 ##expression heatmap
 #library(RColorBrewer)
 heat_colors <- brewer.pal(6, "YlOrRd")
 #display.brewer.all()
 ?pheatmap
 #run heatmap
 pheatmap( subset_sig_norm, 
          color = heat_colors, 
          cluster_rows = T, 
          show_rownames = T,
          annotation = select(meta, c("outcome", "Type", "wagner_grade", "Batch")),
          #border_color = NA, 
          #fontsize = 10, 
          scale = "row", 
          #fontsize_row = 10, 
          #labels_row = label,
          height = 30)
 

 #box
 
subset_sig_box = sig_norm_counts[c("ENSG00000186832"), ]
subset_sig_box = data.frame(subset_sig_box)%>%
                rownames_to_column(var = "SampleID")%>%
                rename(Expression = subset_sig_box) %>% #rename subset_sig_box to Expression
                mutate(SampleID = str_replace_all(SampleID, "-", ".")) %>% 
                left_join (select(meta_wt,SampleID,Type,outcome), by = "SampleID") #only add type and outcome to the merged dataframe
            


ggplot(subset_sig_box, aes(x =outcome, y = Expression, fill = outcome)) +
  geom_boxplot() +
  labs(title = "Gene - KRT16:Keratin16",
       x = "Grouping Variable",
       y = "Expression level") +
  theme_minimal()

 
 
 
 
 


##  Transformation function
vsd <- varianceStabilizingTransformation(dds_2, blind = FALSE)

## Principal component plot of the samples
plotPCA(vsd, intgroup="outcome")

#remove batch effect
mat <- assay(vsdata)
mm <- model.matrix(~outcome, colData(vsdata))
mat <- limma::removeBatchEffect(mat, batch=vsdata$Batch, design=mm)
assay(vsdata) <- mat
plotPCA(vsdata, intgroup = "Batch", pcsToUse = 3:4)














# Top 100 genes based on mean expression
# Calculate the mean expression for each gene
mean_expression <- rowMeans(counts(dds))

# Order genes by mean expression in descending order
ordered_genes <- names(sort(mean_expression, decreasing = TRUE))

# Get the top 100 genes
top_100_genes <- ordered_genes[1:100]

# Display or use the list of top 100 genes
print(top_100_genes)

top = as.data.frame(top_100_genes)




#PCA play
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('PCAtools')

library(PCAtools)


vst <- assay(vsdata)
p <- pca(vst, metadata = meta, removeVar = 0.1)

screeplot(p, axisLabSize = 18, titleLabSize = 22)

biplot(p)
biplot(p, showLoadings = TRUE, lab = NULL)


#Colour by a metadata factor, use a custom label, add lines through origin, and add legend
biplot(p,
       lab = paste0(p$metadata$outcome),
       colby = 'Patient',
       hline = 0, vline = 0,
       legendPosition = 'right')

biplot(p,
       colby = 'outcome', colkey = c('Amputated' = 'forestgreen', 'Healed' = 'purple'),
       # ellipse config
       ellipse = TRUE,
       ellipseType = 't',
       ellipseLevel = 0.95,
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 1.0,
       xlim = c(-125,125), ylim = c(-50, 80),
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)


