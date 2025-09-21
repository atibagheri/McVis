suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(openxlsx)
  library(readr)
  library(dplyr)
  library(zip)
  library(base64enc)
})

# ---------- Helpers ----------
get_org <- function(species) {
  if (identical(species, "mouse")) {
    list(orgdb = org.Mm.eg.db, code = "mmu")
  } else if (identical(species, "human")) {
    list(orgdb = org.Hs.eg.db, code = "hsa")
  } else {
    stop("Unsupported species: must be 'mouse' or 'human'")
  }
}

ensure_nonempty_plot <- function(expr, fallback_title) {
  p <- NULL
  try({
    p <- expr
  }, silent = TRUE)
  if (is.null(p)) {
    p <- ggplot() + geom_blank() +
      ggtitle(fallback_title) +
      theme_bw(base_size = 14)
  }
  p
}

save_plot_png <- function(plot, out_path, width = 10, height = 7, dpi = 300) {
  ggsave(out_path, plot = plot, width = width, height = height, dpi = dpi)
  out_path
}

encode_file_b64 <- function(path) {
  base64encode(readBin(path, "raw", n = file.info(path)$size))
}

safe_enrich_barplot <- function(ego, show = 10, title = "Top GO Terms") {
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    return(ggplot() + geom_blank() + ggtitle("No enriched terms") + theme_bw(base_size = 14))
  }
  barplot(ego, showCategory = show) + ggtitle(title) + theme_bw(base_size = 14)
}

safe_cnetplot <- function(ego, show = 5, fc = NULL, title = "Category–Gene Network") {
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    return(ggplot() + geom_blank() + ggtitle("No network (no enriched terms)") + theme_bw(base_size = 14))
  }
  p <- cnetplot(ego, showCategory = show, foldChange = fc) +
    theme_bw(base_size = 14) +
    theme(panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    ggtitle(title)
  p
}

# ========== MODE 1: gene list only ==========
run_go_from_gene_list <- function(input_file, species) {
  message("[gomap] mode=genelist | species=", species, " | file=", input_file)
  org <- get_org(species)

  # input: one SYMBOL per line
  gene_symbols <- read_lines(input_file) |> unique() |> na.omit()
  gene_symbols <- gene_symbols[nchar(gene_symbols) > 0]
  if (length(gene_symbols) == 0) stop("Input file has no gene symbols.")

  gene_df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org$orgdb)
  unmapped <- setdiff(gene_symbols, gene_df$SYMBOL)

  ego <- tryCatch({
    enrichGO(gene = gene_df$ENTREZID,
             ont = "BP",
             OrgDb = org$orgdb,
             pvalueCutoff = 0.05,
             readable = TRUE)
  }, error = function(e) NULL)

  enr_df <- if (!is.null(ego)) as.data.frame(ego@result) else data.frame()

  # --- temp dir & files (auto-clean later) ---
  tdir <- tempfile("gomap_"); dir.create(tdir)
  on.exit(unlink(tdir, recursive = TRUE, force = TRUE), add = TRUE)

  csv_path <- file.path(tdir, "GO_enrichment_results.csv")
  write.csv(enr_df, csv_path, row.names = FALSE)


# Unmapped genes (TXT)
unmapped_path <- file.path(tdir, "unmapped_genes.txt")
writeLines(unmapped, unmapped_path)



  # plots
  barplot_path <- file.path(tdir, "GO_barplot.png")
  bp <- safe_enrich_barplot(ego, 10, "Top 10 GO Terms")
  save_plot_png(bp, barplot_path)

  # zip everything in temp dir into a zip file (still inside temp)
  # --- zip (flat entries at root of the archive) ---
  zip_path <- file.path(tdir, "gomap_results.zip")
  files_to_zip <- c(csv_path, unmapped_path, barplot_path)

  # Prefer zipr with root=tdir so entries are just basenames
  if ("root" %in% names(formals(zip::zipr))) {
    invisible(zip::zipr(
      zipfile = zip_path,
      files   = basename(files_to_zip),
      root    = tdir
    ))
  } else {
    # Fallback for older zip packages: temporarily setwd into tdir
    old <- setwd(tdir); on.exit(setwd(old), add = TRUE)
    invisible(zip::zip(
      zipfile = zip_path,
      files   = basename(files_to_zip)
    ))
  }


  # base64 payloads for UI
  barplot_b64 <- encode_file_b64(barplot_path)
  zip_b64 <- encode_file_b64(zip_path)
  csv_b64      <- encode_file_b64(csv_path)
  unmapped_b64 <- encode_file_b64(unmapped_path)

  # Return and let caller decide to persist anything. We still remove ONLY the temp content when R cleans up the session.
  list(
    n_input_genes = length(gene_symbols),
    n_mapped = nrow(gene_df),
    n_terms = nrow(enr_df),
    barplot_base64 = barplot_b64,
    cnetplot_base64 = NULL,
    csv_base64      = csv_b64,
    csv_filename    = "GO_enrichment_results.csv",
    unmapped_base64 = unmapped_b64,
    unmapped_filename = "unmapped_genes.txt",
    zip_base64 = zip_b64,
    zip_filename = "gomap_genelist_results.zip"
  )
}

# ========== MODE 2: gene list + log2FC ==========
run_go_from_gene_list_fc <- function(input_file, species) {
  message("[gomap] mode=genelist_fc | species=", species, " | file=", input_file)
  org <- get_org(species)

  deg_data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("Gene.ID", "log2FC") %in% colnames(deg_data))) {
    stop("Input TXT must have columns: Gene.ID and log2FC")
  }
  deg_data <- deg_data |> distinct() |> na.omit()

  gene_df <- bitr(deg_data$Gene.ID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org$orgdb)
  unmapped <- setdiff(deg_data$Gene.ID, gene_df$SYMBOL)
  deg_merged <- merge(gene_df, deg_data, by.x = "SYMBOL", by.y = "Gene.ID")

  if (nrow(deg_merged) == 0) stop("No genes could be mapped to ENTREZ IDs.")

  fc_vector <- deg_merged$log2FC
  names(fc_vector) <- deg_merged$ENTREZID

  ego <- tryCatch({
    enrichGO(gene = deg_merged$ENTREZID,
             OrgDb = org$orgdb,
             keyType = "ENTREZID",
             ont = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             readable = TRUE)
  }, error = function(e) NULL)

  enr_df <- if (!is.null(ego)) as.data.frame(ego@result) else data.frame()

  # --- temp dir & files ---
  tdir <- tempfile("gomap_"); dir.create(tdir)

  csv_path <- file.path(tdir, "GO_enrichment_results.csv")
  write.csv(enr_df, csv_path, row.names = FALSE)

  # Unmapped genes (TXT)
  unmapped_path <- file.path(tdir, "unmapped_genes.txt")
  writeLines(unmapped, unmapped_path)


  # plots
  barplot_path <- file.path(tdir, "GO_barplot.png")
  bp <- safe_enrich_barplot(ego, 10, "Top 10 GO Terms")
  save_plot_png(bp, barplot_path)

  cnet_path <- file.path(tdir, "GO_cnetplot_top5.png")
  cp <- safe_cnetplot(ego, 5, fc_vector, "Top 5 GO terms — gene network")
  save_plot_png(cp, cnet_path, width = 12, height = 9, dpi = 300)

  # --- zip (flat entries at root of the archive) ---
  zip_path <- file.path(tdir, "gomap_results.zip")
  files_to_zip <- c(csv_path, unmapped_path, barplot_path, cnet_path)

  if ("root" %in% names(formals(zip::zipr))) {
    invisible(zip::zipr(
      zipfile = zip_path,
      files   = basename(files_to_zip),
      root    = tdir
    ))
  } else {
    old <- setwd(tdir); on.exit(setwd(old), add = TRUE)
    invisible(zip::zip(
      zipfile = zip_path,
      files   = basename(files_to_zip)
    ))
  }

  # base64 payloads
  barplot_b64  <- encode_file_b64(barplot_path)
  cnet_b64     <- encode_file_b64(cnet_path)
  zip_b64      <- encode_file_b64(zip_path)
  csv_b64      <- encode_file_b64(csv_path)
  unmapped_b64 <- encode_file_b64(unmapped_path)

  list(
    n_input_genes = nrow(deg_data),
    n_mapped = nrow(deg_merged),
    n_terms = nrow(enr_df),
    barplot_base64  = barplot_b64,
    cnetplot_base64 = cnet_b64,
    csv_base64      = csv_b64,
    csv_filename    = "GO_enrichment_results.csv",
    unmapped_base64 = unmapped_b64,
    unmapped_filename = "unmapped_genes.txt",          
    zip_base64      = zip_b64,
    zip_filename    = "gomap_genelist_fc_results.zip"
  )
}
