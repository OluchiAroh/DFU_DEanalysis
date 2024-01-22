if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DOSE")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")

install.packages("ggplot2")


# Load libraries
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

#annotate all genes and module genes

library(annotables)
library(dplyr)
grch38
grch = grch38[, c("ensgene", "symbol", "description")] 
grch_f = distinct(grch)
nrow(grch38)
nrow(grch_f)

#transpose norm.counts
norm.counts_df = t(norm.counts)

#annotate all the genes used for WGCNA 
all_res <- data.frame(norm.counts_df) %>%
  rownames_to_column(var="ensgene") %>% 
  left_join(
    y = grch_f[, c("ensgene", "symbol", "description")],
    by = "ensgene")
nrow(all_res)


#subset the genes present in midnight blue module
mod_genes = gene_module_midnight$gene
all_res_sig <- all_res[all_res$ensgene %in% mod_genes, ]


##START ANALYSIS
###We will be using clusterProfiler to perform over-representation analysis on GO terms associated with our list of significant genes.
#All genes
nrow(all_res) 

#significant genes
nrow(all_res_sig)

#create background data for the test using all genes
all_genes <- as.character(all_res$ensgene)
nrow(all_genes)

#extract significant genes
sig_genes <- as.character(all_res_sig$ensgene)
nrow(sig_genes)

#Run GO enrichment - with all genes in my result

Go_WGCNA_MF <- enrichGO(gene = sig_genes, 
                   universe = all_genes,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "MF", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

#view result
cluster_MF_summary <- data.frame(Go_WGCNA_MF)
write.csv(cluster_go_summary, "results/cluster_summary_tissue_WCGNA")

nrow(as.data.frame(Go_WGCNA_MF))

#remove.packages("yulab.utils")
#install.packages("yulab.utils", dependencies = TRUE)

#plot
D = plot(barplot(Go_WGCNA_MF, showCategory = 25))
D + ggtitle('Enriched GO Molecular Function - Module MidnightBlue:Tissue') 

dotplot(Go_WGCNA_MF, showCategory=25, font.size = 10)

?dotplot


#Group similar terms together
all_ego = enrichplot::pairwise_termsim(Go_WGCNA_BP)

## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(all_ego, showCategory = 70, font.size = 10)



#GSEA - utilizes the gene-level statistics or log2 fold changes for all genes to look to see whether gene sets for particular biological pathways are enriched among the large positive or negative fold changes.

#GO pathways

#extract column of interest - logfoldchange - 

#order by logfoldchange or stat or any metrics you want

G_res = all_res[order(-all_res$log2FoldChange),]
G_res

gene_list = G_res$log2FoldChange
names(gene_list) = G_res$ensgene

?gseGO
gse <- gseGO(gene_list,
             ont = "BP",
             keyType = "ENSEMBL",
             OrgDb = "org.Hs.eg.db",
             pvalueCutoff = 0.05,
             eps = 1e-300)


as.data.frame(gse)
Go_path <- data.frame(gse)
write.csv(Go_path, "results/Go_path_swab")

gseaplot(gse, geneSetID = 1)

#KEGG PATHWAY
grch_a = grch38[, c("entrez","ensgene")] 
grch_tb = distinct(grch_a)

nrow(grch38)
nrow(grch_tb)

all_res_en <- data.frame(resLFC) %>%
  rownames_to_column(var="ensgene") %>% 
  left_join(
    y = grch_tb[, c("entrez", "ensgene")],
    by = "ensgene")

nrow(all_res_en)

#remove NAs
## Remove any NA values (reduces the data by quite a bit)
res_entrez <- dplyr::filter(all_res_en, entrez != "NA")
nrow((res_entrez))

## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrez) == F), ]

## Extract the foldchanges
foldchanges <- res_entrez$log2FoldChange

## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$entrez

## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)
head(foldchanges)

set.seed(123456)

## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    minGSSize = 5, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = T,
                    eps = 1e-300)


gseaKEGG_results <- gseaKEGG@result
KEGG_path <- data.frame(gseaKEGG)
write.csv(KEGG_path, "results/KEGG_path_swab")

## Plot the GSEA plot for a single enriched pathway, `hsa03040`
gseaplot(gseaKEGG, geneSetID = 'hsa04740')


#pathway view

detach("package:dplyr", unload=TRUE) # first unload dplyr to avoid conflicts

## Output images for a single significant KEGG pathway
pathview(gene.data = foldchanges,
         pathway.id = "hsa04740",
         species = "hsa",
         limit = list(gene = 2, # value gives the max/min limit for foldchanges
                      cpd = 1))


## Output images for all significant KEGG pathways
#get_kegg_plots <- function(x) {
#  pathview(gene.data = foldchanges, 
#          pathway.id = gseaKEGG_results$ID[x], 
#           species = "hsa",
#           limit = list(gene = 2, cpd = 1))
#}

#purrr::map(1:length(gseaKEGG_results$ID), 
#           get_kegg_plots)



#category netplot -shows genes associated with the top five most significant GO TERM

fold_change = all_res_sig$log2FoldChange

names(fold_change) <- all_res_sig$ensgene

#Cnetplot
q = cnetplot(Go_WGCNA_BP, 
             categorySize="pvalue", 
             showCategory = 25, 
             foldChange=fold_change, 
             vertex.label.font=6)

?cnetplot
q$labels


rownames(all_res_sig$ensgene)
all_res_sig$ensgene
