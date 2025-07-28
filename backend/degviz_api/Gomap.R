library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(pathview)
library(ggplot2)
library(openxlsx)
library(readr)
library(dplyr)
library(zip)
library(base64enc)

# Helper: pick OrgDb and organism code
get_org <- function(species) {
  if (species == "mouse") {
    list(orgdb = org.Mm.eg.db, GO = "mmu")
  } else if (species == "human") {
    list(orgdb = org.Hs.eg.db, GO = "hsa")
  } else {
    stop("Unsupported species: must be 'mouse' or 'human'")
  }
}

##############################################
# MODE 1: Only gene list
##############################################
run_go_from_gene_list <- function(input_file, species) {
  message("[run_gomap] Starting analysis...")
  message("📄 input_file: ", input_file)
  message("🔧 mode: genelist")
  message("🧬 species: ", species)
  
  output_dir <- "./gomap_output"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  org <- get_org(species)
  
  # Read gene symbols
  gene_symbols <- read_lines(input_file) %>% unique() %>% na.omit()
  message("✅ Columns: Gene.ID")
  
  gene_df <- bitr(gene_symbols,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org$orgdb)
  message("✅ Mapped genes: ", nrow(gene_df))
  
  # Save unmapped
  unmapped <- setdiff(gene_symbols, gene_df$SYMBOL)
  writeLines(unmapped, file.path(output_dir, "unmapped_genes.txt"))
  
  # Enrichment
  go_enrich <- enrichGO(
    gene = gene_df$ENTREZID,
    ont = "BP",
    OrgDb = org$orgdb,
    pvalueCutoff = 0.05
  )
  
  enrich_df <- as.data.frame(go_enrich@result)
  message("✅ Enrichment results: ", nrow(enrich_df))
  csv_path <- file.path(output_dir, "GO_enrichment_results.csv")
  write.csv(enrich_df, file = csv_path, row.names = FALSE)
  
  # Barplot
  barplot_file <- file.path(output_dir, "GO_barplot.png")
  bp <- barplot(go_enrich, showCategory = 10) +
    ggtitle("Top 10 GO Terms") +
    theme_bw(base_size = 14)
  ggsave(barplot_file, plot = bp, width = 10, height = 7, dpi = 300)
  message("✅ Saved barplot: ", barplot_file)
  
  # Base64 for UI
  barplot_base64 <- base64encode(readBin(barplot_file, "raw", file.info(barplot_file)$size))
  
  # Zip results
  zip_file <- file.path(output_dir, "gomap_results.zip")
  if (file.exists(zip_file)) file.remove(zip_file)
  zip(zip_file, files = list.files(output_dir, full.names = TRUE))
  message("📦 Zip created: ", zip_file)
  
  return(list(
    barplot_base64 = barplot_base64,
    cnetplot_base64 = NULL,
    unmapped = file.path(output_dir, "unmapped_genes.txt"),
    csv_path = file.path(output_dir, "GO_enrichment_results.csv"),
    zip_file = basename(zip_file),
    output_dir = output_dir
  ))
}


##############################################
# MODE 2: Gene list + fold change
##############################################
run_go_from_gene_list_fc <- function(input_file, species) {
  message("[run_gomap] Starting analysis...")
  message("📄 input_file: ", input_file)
  message("🔧 mode: genelist_fc")
  message("🧬 species: ", species)
  
  output_dir <- "./gomap_output"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  org <- get_org(species)
  
  # Read FC data
  deg_data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (!all(c("Gene.ID", "log2FC") %in% colnames(deg_data))) {
    stop("❌ Input TXT must have columns: Gene.ID and log2FC")
  }
  deg_data <- deg_data %>% na.omit() %>% unique()
  
  gene_df <- bitr(deg_data$Gene.ID,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org$orgdb)
  unmapped_genes <- setdiff(deg_data$Gene.ID, gene_df$SYMBOL)
  writeLines(unmapped_genes, file.path(output_dir, "unmapped_genes.txt"))
  
  deg_merged <- merge(gene_df, deg_data, by.x = "SYMBOL", by.y = "Gene.ID")
  message("✅ Mapped genes: ", nrow(deg_merged))
  
  fc_vector <- deg_merged$log2FC
  names(fc_vector) <- deg_merged$ENTREZID
  message("✅ Fold-change vector: ", length(fc_vector))
  
  go_enrich <- enrichGO(
    gene = deg_merged$ENTREZID,
    OrgDb = org$orgdb,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  enrich_df <- as.data.frame(go_enrich@result)
  csv_path <- file.path(output_dir, "GO_enrichment_results.csv")
  write.csv(enrich_df, file = csv_path, row.names = FALSE)
  
  # Barplot
  barplot_file <- file.path(output_dir, "GO_barplot.png")
  bp <- barplot(go_enrich, showCategory = 10) +
    ggtitle("Top 10 GO Terms") +
    theme_bw(base_size = 14)
  ggsave(barplot_file, plot = bp, width = 10, height = 7, dpi = 300)
  message("✅ Saved barplot: ", barplot_file)
  
  barplot_base64 <- base64encode(readBin(barplot_file, "raw", file.info(barplot_file)$size))
  
  # Cnetplot
  cnetplot_file <- file.path(output_dir, "GO_cnetplot_top5.png")
  cp <- cnetplot(go_enrich, showCategory = 5, foldChange = fc_vector) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  ggsave(cnetplot_file, plot = cp, width = 12, height = 9, dpi = 300)
  message("✅ Saved cnetplot: ", cnetplot_file)
  
  cnetplot_base64 <- base64encode(readBin(cnetplot_file, "raw", file.info(cnetplot_file)$size))
  
  
   # zip
    zip_file <- file.path(output_dir, "gomap_results.zip")
    if (file.exists(zip_file)) file.remove(zip_file)
    zip::zip(zipfile = zip_file, files = list.files(output_dir, full.names = TRUE))

    return(list(
      barplot_base64 = barplot_base64,
      cnetplot_base64 = cnetplot_base64,
      unmapped = file.path(output_dir, "unmapped_genes.txt"), 
      csv_path = file.path(output_dir, "GO_enrichment_results.csv"),
      zip_file = basename(zip_file),
      output_dir = normalizePath(output_dir)
  
    ))
}
