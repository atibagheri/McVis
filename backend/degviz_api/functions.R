# ============================================
# ✅ functions.R
# Contains both Venn/UpSet and Enrichment functions
# ============================================
library(base64enc)

generate_gene_overlap_plot <- function(file_paths, file_labels = NULL, output_path = "gene_overlap_plot.png") {
  if (length(file_paths) < 2 || length(file_paths) > 7) {
    stop("Please provide between 2 and 7 gene list files.")
  }

  if (is.null(file_labels)) {
    file_labels <- tools::file_path_sans_ext(basename(file_paths))
  }

  gene_lists <- lapply(file_paths, function(f) {
    df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (!"Gene.ID" %in% colnames(df)) stop(paste("Missing 'Gene.ID' in", f))
    df$Gene.ID
  })
  names(gene_lists) <- file_labels

  if (length(gene_lists) <= 5) {
    # Venn
    library(VennDiagram)
    futile.logger::flog.threshold(futile.logger::ERROR)
    venn <- venn.diagram(
      x = gene_lists,
      filename = NULL,
      imagetype = "png",
      fill = c("red","blue","green","orange","purple","yellow","pink")[1:length(gene_lists)],
      alpha = 0.5,
      cat.cex = 1.2,
      margin = 0.1
    )
    png(output_path, width = 1200, height = 1000, res = 200)
    grid::grid.draw(venn)
    dev.off()
  } else {
    # UpSet
    library(UpSetR)
    all_genes <- unique(unlist(gene_lists))
    gene_matrix <- sapply(gene_lists, function(g) all_genes %in% g)
    gene_df <- as.data.frame(gene_matrix)
    gene_df$Gene <- all_genes
    colnames(gene_df) <- c(file_labels, "Gene")
    png(output_path, width = 1600, height = 1200, res = 200)
    upset(gene_df, sets = file_labels, order.by = "freq")
    dev.off()
  }

  message("✅ Overlap plot saved to: ", output_path)
  # After saving the plot, read the image as raw bytes and return as base64 string
  img_bytes <- readBin(output_path, what = "raw", n = file.info(output_path)$size)
  img_base64 <- base64enc::base64encode(img_bytes)
  return(list(
    output_path = output_path,
    png_base64 = img_base64
  ))
}


#Enrichment functions
# ============================================

run_enrichment_analysis <- function(filePath, species = "mouse") {
  print(getwd())

  if (!file.exists(filePath)) stop("❌ Input file does not exist: ", filePath)

  species <- tolower(species)
  if (!(species %in% c("human", "mouse"))) {
    stop("❌ Invalid species: ", species, ". Use 'human' or 'mouse'.")
  }

  gene_list <- read.table(filePath, header = TRUE, stringsAsFactors = FALSE)
  if (!"Gene.ID" %in% colnames(gene_list)) stop("❌ Column 'Gene.ID' not found.")

  genes <- as.character(gene_list$Gene.ID)

  if (species == "human") {
    OrgDb <- org.Hs.eg.db
    kegg_org <- "hsa"
  } else {
    OrgDb <- org.Mm.eg.db
    kegg_org <- "mmu"
  }

  entrez <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb)
  if (is.null(entrez) || nrow(entrez) == 0) stop("❌ No valid genes mapped.")

  go_res <- enrichGO(gene = entrez$ENTREZID, OrgDb = OrgDb, ont = "BP", pvalueCutoff = 0.05)
  kegg_res <- enrichKEGG(gene = entrez$ENTREZID, organism = kegg_org, pvalueCutoff = 0.05)

  output_dir <- "./Enrichment"
  go_pdf <- file.path(output_dir, "go_enrichment_plot.pdf")
  kegg_pdf <- file.path(output_dir, "kegg_enrichment_plot.pdf")
  go_txt <- file.path(output_dir, "go_enrichment_results.txt")
  kegg_txt <- file.path(output_dir, "kegg_enrichment_results.txt")
  # Save plots
  pdf(go_pdf, width = 14, height = 10); print(dotplot(go_res, showCategory = 15)); dev.off()
  pdf(kegg_pdf, width = 14, height = 10); print(dotplot(kegg_res, showCategory = 15)); dev.off()

  # Save enrichment tables as TXT (tab‑separated)
  write.table(as.data.frame(go_res),
              file = go_txt,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  write.table(as.data.frame(kegg_res),
              file = kegg_txt,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)



  message("✅ Enrichment completed for species: ", species)
  message("📂 Results saved in: ", normalizePath(output_dir))

  # Return all file paths + directory
  return(list(
    output_dir = normalizePath(output_dir),
    go_pdf = normalizePath(go_pdf),
    kegg_pdf = normalizePath(kegg_pdf),
    go_txt = normalizePath(go_txt),
    kegg_txt = normalizePath(kegg_txt)
  ))
}
