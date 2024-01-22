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

# Load libraries
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

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
nrows(sig_genes)

#Run GO enrichment - with all genes in my result

Go_res <- enrichGO(gene = sig_genes, 
                universe = all_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

#view result
cluster_go_summary <- data.frame(Go_res)
write.csv(cluster_go_summary, "results/cluster_summary_tissue")

nrow(as.data.frame(Go_res))

#plot
D = plot(barplot(Go_res, showCategory = 50))
D + ggtitle('Enriched GO Biological Process - All significant genes:Tissue') 

dotplot(Go_res, showCategory=25, font.size = 10)

?dotplot

##Focusing only on upregulated genes OR DOWNREGULATED
sig_up <- dplyr::filter(all_res, padj < 0.05 & log2FoldChange > 0)
sig_down <- dplyr::filter(all_res, padj < 0.05 & log2FoldChange < 0)

#dataframe it
sig_up_genes <- as.character(sig_up$ensgene)
sig_down_genes <- as.character(sig_down$ensgene)

#upregulated
Go_res_up <- enrichGO(gene = sig_up_genes, 
                   universe = all_genes,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "ALL", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

nrow(as.data.frame(Go_res_up))
#plot
C = plot(barplot(Go_res_up, showCategory = 25))
C + ggtitle('Enriched GO Biological process - Upregulated genes:Tissue') 

dotplot(Go_res_up, showCategory=25, font.size = 10)

up_go_summary <- data.frame(Go_res_up)
write.csv(up_go_summary, "results/up_go_summary_tissue")

#downregulated
Go_res_down <- enrichGO(gene = sig_down_genes, 
                   universe = all_genes,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

as.data.frame(Go_res_down)
#view result
down_go_summary <- data.frame(Go_res_down)
write.csv(down_go_summary, "results/down_go_summary_swab")

#plot
E = plot(barplot(Go_res_down, showCategory = 15))
E + ggtitle('Enriched GO Biological process - Downregulated genes:Swab') 

dotplot(Go_res_down, showCategory=50)

#Group similar terms together
all_ego = enrichplot::pairwise_termsim(Go_res)
up_ego <- enrichplot::pairwise_termsim(Go_res_up)

## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(all_ego, showCategory = 50, font.size = 10)
emapplot(up_ego, showCategory = 50, font.size = 10)


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
q = cnetplot(up_ego, 
        categorySize="pvalue", 
        showCategory = 25, 
        foldChange=fold_change, 
        vertex.label.font=6)

?cnetplot
q$labels

## If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
#OE_foldchanges <- ifelse(OE_foldchanges > 2, 2, OE_foldchanges)
#OE_foldchanges <- ifelse(OE_foldchanges < -2, -2, OE_foldchanges)

#cnetplot(ego, 
#        categorySize="pvalue", 
#         showCategory = 5, 
#       foldChange=OE_foldchanges, 
#        vertex.label.font=6)

## To show significant processess not among the top five
#Subset the "ego" variable without overwriting original `ego` variable - 
#ego2 <- ego

#ego2@result <- ego@result[c(1,3,4,8,9),]

## Plotting terms of interest
#cnetplot(ego2, 
#         categorySize="pvalue", 
#         foldChange=OE_foldchanges, 
#         showCategory = 5, 
 #        vertex.label.font=6)
