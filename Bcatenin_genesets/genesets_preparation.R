# Preparing genesets

# Loading libraries ----
library(tidyverse)
library(dplyr) # used for data manipulation
library(biomaRt) # for gene id conversion functions
library(RCurl) # if proxy issues

library(systemPipeR) # used to generate vennplots
library(UpSetR) # used to visualize overlaps between datasets

# Setting main directories ----
setwd("./Bcatenin_genesets/")

DATA_DIR="./data/"
GENESETS_DIR="./genesets/"

# Defining functions ----
# for conversion of human to mouse ids
# functional enrichment
queryBiomartMouse2human <- function(geneList = NULL,
                                    filters = "ensembl_gene_id",
                                    ensembl_version = ensembl_version,
                                    biomart_host = biomart_host) {
  # Query biomart for mouse ensembl_ids and return human orthologues and attributes
  
  human <- biomaRt::useMart("ensembl",
                            dataset = "hsapiens_gene_ensembl",
                            host = biomart_host,
                            version = ensembl_version)
  mouse <- biomaRt::useMart("ensembl",
                            dataset = "mmusculus_gene_ensembl",
                            host = biomart_host,
                            version = ensembl_version)
  out <- biomaRt::getLDS(
    attributes = c("ensembl_gene_id"),
    filters = filters,
    values = geneList,
    mart = mouse,
    attributesL = c(
      "ensembl_gene_id",
      "hgnc_symbol",
      "entrezgene",
      "description",
      "chromosome_name",
      "start_position",
      "end_position",
      "strand"
    ),
    martL = human
  )
  return(out)
}

convertIDs_mouse2human <- function(geneList = NULL,
                                   ensembl_version = "Ensembl Genes 92",
                                   biomart_host = "http://apr2018.archive.ensembl.org") {
  require("dplyr")
  require("biomaRt")
  require("RCurl")
  # use proxy if needed
  options(RCurlOptions = list(proxy = "wwwcache.gla.ac.uk:8080",
                              http.version = HTTP_VERSION_1_0))
  
  cat("Converting mouse ensembl_ids to human ensembl_ids.\n")
  
  ensembl_version = ensembl_version #"Ensembl Genes 92"
  cat("Using ensembl version:", ensembl_version, "\n")
  
  # Converting mouse to human entrez_ids
  human_ids <- queryBiomartMouse2human(
      geneList = geneList,
      ensembl_version = ensembl_version,
      biomart_host = biomart_host
    )
  
  # Renaming default column names
  human_ids <- human_ids %>%
    dplyr::rename(
      ensembl_id_mouse = Gene.stable.ID,
      ensembl_id = Gene.stable.ID.1,
      hgnc_symbol = HGNC.symbol,
      hgnc_description = Gene.description,
      entrezgene = NCBI.gene.ID,
      chromosome_name = Chromosome.scaffold.name,
      start_position = Gene.start..bp.,
      end_position = Gene.end..bp.,
      strand = Strand
    ) %>%
    dplyr::mutate(entrezgene = as.character(entrezgene)) %>%
    as.data.frame(.)

  # printing basic summary
  cat("Mouse ensembl_ids:", length(geneList), "\n")
  cat("  converted to human: \n")
  cat("  ensembl_ids:", length(human_ids$ensembl_id), "( duplicates:", sum(duplicated(human_ids$ensembl_id)), ")\n")
  cat("  entrez_ids:", length(human_ids$entrezgene), " ( duplicates:", sum(duplicated(human_ids$entrezgene)), ")\n")
  cat("  hgnc_symbols:", length(human_ids$hgnc_symbol), "( duplicates:", sum(duplicated(human_ids$hgnc_symbol)), ")\n")

  return(human_ids)
}

## Extracting and comparing significantly differentially expressed genes between conditions ----
# gemm_orthografts_signif_results.RData:
#  file contains significantly DE genes for orthografts and GEMM as lists of UP and DOWN regulated genes for each of the comparisons
#  it also contains ensembl_ids of all detected genes (across conditions) used as a background for enrichment analysis
# genes were filtered based on (padj_cutoff=0.05, log2FC_cutoff=1) using the following:
# orthografts_CP1vsP1_signif <-  orthografts_CP1vsP1 %>%
#   filter((!is.na(padj) & padj < padj_cutoff) & abs(log2FoldChange) > log2FC_cutoff) 
load(file = paste0(DATA_DIR, "gemm_orthografts_signif_results.RData"))

venn_GEMM_orthografts_list = list(CP1vsP1_UP = as.character(orthografts_CP1vsP1_signif$UP),
                                        CP2vsP1_UP = as.character(orthografts_CP2vsP1_signif$UP),
                                        CPhetvsPhom_UP = as.character(GEMM_CPhetvsPhom_signif$UP),
                                        CP1vsP1_DOWN = as.character(orthografts_CP1vsP1_signif$DOWN),
                                        CP2vsP1_DOWN = as.character(orthografts_CP2vsP1_signif$DOWN),
                                        CPhetvsPhom_DOWN = as.character(GEMM_CPhetvsPhom_signif$DOWN))

# Generating overlaps between datasets ----
interset_GEMM_orthografts <- systemPipeR::overLapper(venn_GEMM_orthografts_list, type="intersects")

interset_GEMM_orthografts_df <- as.data.frame(intersectmatrix(interset_GEMM_orthografts)) %>%
  tibble::rownames_to_column(var="ensembl_id") %>%
  dplyr::mutate(all_UP=(CP1vsP1_UP+CP2vsP1_UP+CPhetvsPhom_UP),
         all_DOWN=(CP1vsP1_DOWN+CP2vsP1_DOWN+CPhetvsPhom_DOWN)) # summing UP and DOWN

interset_GEMM_orthografts_df_UP <- interset_GEMM_orthografts_df %>%
  filter(all_UP == 3)
interset_GEMM_orthografts_df_DOWN <- interset_GEMM_orthografts_df %>%
  filter(all_DOWN == 3)

# overlap across all 3 datasets (2 orthografts and GEMM)
interset_GEMM_orthografts_WNT_list <- venn_GEMM_orthografts_list
interset_GEMM_orthografts_WNT_list$beta_catenin_UP <- as.character(interset_GEMM_orthografts_df_UP$ensembl_id)
interset_GEMM_orthografts_WNT_list$beta_catenin_DOWN <- as.character(interset_GEMM_orthografts_df_DOWN$ensembl_id)

interset_GEMM_orthografts_WNT <- systemPipeR::overLapper(interset_GEMM_orthografts_WNT_list, type="intersects")
interset_GEMM_orthografts_WNT_df <- as.data.frame(systemPipeR::intersectmatrix(interset_GEMM_orthografts_WNT))

# visualize overlaps
# library(UpSetR)
# fix naming with '-' otherwise upset crashes
colnames(interset_GEMM_orthografts_WNT_df) <- gsub(pattern = "-", replacement = "_", colnames(interset_GEMM_orthografts_WNT_df))
#colnames(interset_GEMM_orthografts_WNT_df) <- gsub(pattern = "_", replacement = "", colnames(interset_GEMM_orthografts_WNT_df))
colnames(interset_GEMM_orthografts_WNT_df) <- gsub(pattern = " ", replacement = "_", colnames(interset_GEMM_orthografts_WNT_df))

# reorder sets
sets_order <- rev(c("beta_catenin_UP", "CPhetvsPhom_UP", "CP1vsP1_UP", "CP2vsP1_UP",
                    "beta_catenin_DOWN", "CPhetvsPhom_DOWN", "CP1vsP1_DOWN", "CP2vsP1_DOWN"))
interset_GEMM_orthografts_WNT_df <- interset_GEMM_orthografts_WNT_df[, match(sets_order, colnames(interset_GEMM_orthografts_WNT_df))]


# generating figure and saving as pdf
pdf(file = paste0(GENESETS_DIR,"WNT_ensemblID_genesets_intersections.pdf"), onefile=FALSE)
UpSetR::upset(interset_GEMM_orthografts_WNT_df, 
              sets.bar.color = c(rep("#56B4E9", 4), rep("red", 4)),
              queries = list(list(query = intersects, 
                                  params = list("beta_catenin_UP", "CPhetvsPhom_UP", "CP1vsP1_UP", "CP2vsP1_UP"), active = T),
                             list(query = intersects, 
                                  params = list("beta_catenin_DOWN", "CPhetvsPhom_DOWN", "CP1vsP1_DOWN", "CP2vsP1_DOWN"), active = T)),
              keep.order=TRUE, 
              sets = colnames(interset_GEMM_orthografts_WNT_df),
              order.by = c("degree", "freq"), decreasing = c(F, T),
              text.scale = 0.7) #nintersects = 20
dev.off()

## Generating genesets (gmt) files ----
# Converting mouse ensembl_ids to human hgnc gene symbols ----
# upregulated genes
interset_GEMM_orthografts_df_UP_human <- convertIDs_mouse2human(geneList = interset_GEMM_orthografts_df_UP$ensembl_id, 
                                                                ensembl_version="Ensembl Genes 92", 
                                                                biomart_host="http://apr2018.archive.ensembl.org")

# Converting mouse ensembl_ids to human ensembl_ids.
# Using ensembl version: Ensembl Genes 92 
# Mouse ensembl_ids: 175 
# converted to human: 
#   ensembl_ids: 147 ( duplicates: 4 )
# entrez_ids: 147  ( duplicates: 5 )
# hgnc_symbols: 147 ( duplicates: 4 )

# downregulated genes
interset_GEMM_orthografts_df_DOWN_human <- convertIDs_mouse2human(geneList = interset_GEMM_orthografts_df_DOWN$ensembl_id, 
                                                                  ensembl_version="Ensembl Genes 92", 
                                                                  biomart_host="http://apr2018.archive.ensembl.org")

# Converting mouse ensembl_ids to human ensembl_ids.
# Using ensembl version: Ensembl Genes 92 
# Mouse ensembl_ids: 650 
# converted to human: 
#   ensembl_ids: 660 ( duplicates: 30 )
# entrez_ids: 660  ( duplicates: 37 )
# hgnc_symbols: 660 ( duplicates: 37 )

# cleaning human orthologues
interset_GEMM_orthografts_df_UP_human_symbols <- interset_GEMM_orthografts_df_UP_human %>%
  dplyr::filter(!is.na(hgnc_symbol)) %>%
  dplyr::filter(hgnc_symbol != "") %>%
  dplyr::filter(!duplicated(hgnc_symbol)) %>%
  dplyr::arrange(hgnc_symbol) # order genes alphabetically

interset_GEMM_orthografts_df_DOWN_human_symbols <- interset_GEMM_orthografts_df_DOWN_human %>%
  dplyr::filter(!is.na(hgnc_symbol)) %>%
  dplyr::filter(hgnc_symbol != "") %>%
  dplyr::filter(!duplicated(hgnc_symbol)) %>%
  dplyr::arrange(hgnc_symbol) # order genes alphabetically

GEMM_orthografts_UP_human_symbols_geneset_tsv <- paste(interset_GEMM_orthografts_df_UP_human_symbols$hgnc_symbol, collapse = "\t")
GEMM_orthografts_DOWN_human_symbols_geneset_tsv <- paste(interset_GEMM_orthografts_df_DOWN_human_symbols$hgnc_symbol, collapse = "\t")

# saving extended WTN genesets ----
Bcatenin_genesets <- data.frame(geneset=c("beta-catenin UP",
                                     "beta-catenin DOWN"),
                           description=c("Genes significantly upregulated in CP1_vs_P1, CP2_vs_P1 and CPhet_vs_Phom",
                                         "Genes significantly downregulated in CP1_vs_P1, CP2_vs_P1 and CPhet_vs_Phom"),
                           genes=c(GEMM_orthografts_UP_human_symbols_geneset_tsv,
                                   GEMM_orthografts_DOWN_human_symbols_geneset_tsv))

write.table(Bcatenin_genesets, file=paste0(GENESETS_DIR, "Bcatenin_genesets.gmt"), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)