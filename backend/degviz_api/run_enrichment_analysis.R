#!/usr/bin/env Rscript

# ================================
# Enrichment Analysis Script
# ================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_enrichment_analysis.R <input_file> <output_dir> <species: mouse|human>")
}
input_file <- args[1]
output_dir <- args[2]
species <- tolower(args[3])

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Function to install missing packages
install_if_missing <- function(pkg, bioconductor = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioconductor) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Install required packages
core_pkgs <- c("ggplot2", "stringr")
bio_pkgs <- c("clusterProfiler", "enrichplot", "DOSE","org.Hs.eg.db","org.Mm.eg.db")
for (pkg in core_pkgs) install_if_missing(pkg)
for (pkg in bio_pkgs) install_if_missing(pkg, bioconductor = TRUE)

# Set species-specific database
if (species == "mouse") {
  install_if_missing("org.Mm.eg.db", bioconductor = TRUE)
  library(org.Mm.eg.db)
  OrgDb <- org.Mm.eg.db
  organism_code <- "mmu"
} else if (species == "human") {
  install_if_missing("org.Hs.eg.db", bioconductor = TRUE)
  library(org.Hs.eg.db)
  OrgDb <- org.Hs.eg.db
  organism_code <- "hsa"
} else {
  stop("❌ Species must be either 'mouse' or 'human'")
}

# Read gene list
gene_list <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
genes <- as.character(gene_list$Gene.ID)

# Map to Entrez IDs
entrez_genes <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
unmapped <- setdiff(genes, entrez_genes$SYMBOL)
if (length(unmapped) > 0) {
  cat("⚠️ Unmapped genes:\n")
  print(unmapped)
} else {
  cat("✅ All genes mapped.\n")
}

# ===== GO Enrichment =====
go_enrich <- enrichGO(gene = entrez_genes$ENTREZID,
                      OrgDb = OrgDb,
                      ont = "BP",
                      pvalueCutoff = 0.05,
                      readable = TRUE)
go_enrich@result$Description <- stringr::str_to_title(go_enrich@result$Description)

go_plot <- dotplot(go_enrich, showCategory = 15, x = "Count") +
  ggtitle("GO Enrichment Analysis") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold"))

ggsave(filename = "go_enrichment_plot.pdf", plot = go_plot,
       path = output_dir, width = 14, height = 10)

write.table(as.data.frame(go_enrich),
            file = file.path(output_dir, "go_enrichment_results.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ===== KEGG Enrichment =====
kegg_enrich <- enrichKEGG(gene = entrez_genes$ENTREZID,
                          organism = organism_code,
                          pvalueCutoff = 0.05)
if (nrow(as.data.frame(kegg_enrich)) > 0) {
  kegg_enrich@result$Description <- stringr::str_to_title(kegg_enrich@result$Description)
  
  kegg_plot <- dotplot(kegg_enrich, showCategory = 15, x = "Count") +
    ggtitle("KEGG Enrichment Analysis") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 18, face = "bold"))
  
  ggsave(filename = "kegg_enrichment_plot.pdf", plot = kegg_plot,
         path = output_dir, width = 14, height = 10)
  
  write.table(as.data.frame(kegg_enrich),
              file = file.path(output_dir, "kegg_enrichment_results.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  cat("⚠️ No KEGG enrichment results found.\n")
}

cat("🎉 Analysis complete! Results saved in:", output_dir, "\n")
