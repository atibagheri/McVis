suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(pathview)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(zip)
  library(base64enc)
})

# ---- helpers
get_org <- function(species) {
  if (identical(species, "mouse")) {
    list(orgdb = org.Mm.eg.db, kegg = "mmu")
  } else if (identical(species, "human")) {
    list(orgdb = org.Hs.eg.db, kegg = "hsa")
  } else {
    stop("Unsupported species: 'mouse' or 'human'")
  }
}
b64 <- function(path) base64enc::base64encode(readBin(path, "raw", n = file.info(path)$size))

flat_zip <- function(zip_path, files, root) {
  files <- files[file.exists(files)]
  if (!length(files)) { file.create(zip_path); return(zip_path) }
  if ("root" %in% names(formals(zip::zipr))) {
    invisible(zip::zipr(zipfile = zip_path, files = basename(files), root = root))
  } else {
    old <- setwd(root); on.exit(setwd(old), add = TRUE)
    invisible(zip::zip(zipfile = zip_path, files = basename(files)))
  }
  zip_path
}

# simple patterns
is_base <- function(x, k) {
  # base maps from pathview-native or our _BASE suffix
  grepl(paste0("^", k, "\\d+.*gene.*\\.png$"), basename(x), ignore.case=TRUE) |
    grepl("_BASE.*\\.png$", basename(x), ignore.case=TRUE) |
    grepl("\\.pathview\\.png$", basename(x), ignore.case=TRUE)
}
is_fc <- function(x) {
  grepl("_FC.*\\.png$", basename(x), ignore.case=TRUE) |
    grepl("foldchange.*\\.png$", basename(x), ignore.case=TRUE)
}

# =========================================================
# MODE 1: Gene list only (base highlight maps)
# =========================================================
run_kegg_from_gene_list <- function(input_file, species) {
  org  <- get_org(species)

  # work in a temp sandbox; delete when done
  tdir <- tempfile("kegg_"); dir.create(tdir)
  on.exit(unlink(tdir, recursive = TRUE, force = TRUE), add = TRUE)
  oldwd <- setwd(tdir); on.exit(setwd(oldwd), add = TRUE)

  # read / map
  gene_symbols <- read_lines(input_file) |> unique()
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & nzchar(gene_symbols)]
  if (!length(gene_symbols)) stop("Input file has no gene symbols.")

  gene_df  <- bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org$orgdb)
  unmapped <- setdiff(gene_symbols, gene_df$SYMBOL)
  writeLines(unmapped, file.path(tdir, "unmapped_genes.txt"))

  ek <- tryCatch(enrichKEGG(gene = gene_df$ENTREZID, organism = org$kegg, pvalueCutoff = 0.05),
                 error = function(e) NULL)
  enr_df <- if (!is.null(ek)) as.data.frame(ek@result) else data.frame()
  full_csv <- file.path(tdir, "full_KEGG_enrichment_results.csv")
  write.csv(enr_df, full_csv, row.names = FALSE)

  barplot_path <- NULL
  top_csv      <- NULL

  if (nrow(enr_df) > 0) {
    top10 <- enr_df |> arrange(p.adjust) |> head(10)
    top_csv <- file.path(tdir, "top10_KEGG_enrichment.csv")
    write.csv(top10, top_csv, row.names = FALSE)

    top10$Description <- factor(top10$Description, levels = rev(top10$Description))
    p <- ggplot(top10, aes(x = Count, y = Description, fill = p.adjust)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "#d73027", high = "#4575b4", name = "p.adjust") +
      labs(title = "Top 10 KEGG Pathways", x = "Gene Count", y = NULL) +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            axis.text.y = element_text(size = 12))
    barplot_path <- file.path(tdir, "KEGG_Top10_Barplot.png")
    ggsave(barplot_path, p, width = 10, height = 6, dpi = 300, bg="white")

    # highlight-only maps
    highlight <- setNames(rep(1, length(gene_df$ENTREZID)), gene_df$ENTREZID)
    for (pid in gsub(org$kegg, "", top10$ID)) {
      try(pathview(gene.data  = highlight,
                   pathway.id = pid,
                   species    = org$kegg,
                   kegg.native= TRUE,
                   same.layer = TRUE,
                   gene.idtype= "entrez",
                   limit      = list(gene = 1, cpd = 1),
                   low  = list(gene = NA),
                   mid  = list(gene = NA),
                   high = list(gene = NA),
                   kegg.dir   = tdir), silent = TRUE)
    }
  }

  # collect maps
  all_png <- list.files(tdir, pattern="\\.png$", full.names = TRUE)
  base_maps <- all_png[ is_base(all_png, org$kegg) & !is_fc(all_png) ]

  # zip
  zip_path <- file.path(tdir, "kegg_results.zip")
  files_to_zip <- c(full_csv, if (length(top_csv)) top_csv,
                    file.path(tdir, "unmapped_genes.txt"),
                    if (length(barplot_path)) barplot_path,
                    base_maps)
  flat_zip(zip_path, files_to_zip, root = tdir)

  list(
    success             = TRUE,
    mode                = "genelist",
    species             = species,
    barplot_base64      = if (length(barplot_path)) b64(barplot_path) else NULL,
    maps                = setNames(lapply(base_maps, b64), basename(base_maps)),
    foldchange_maps     = list(),
    zip_base64          = b64(zip_path),
    zip_filename        = "kegg_genelist_results.zip"
  )
}

# =========================================================
# MODE 2: Gene list + log2FC (base + FC maps)
# =========================================================
run_kegg_from_gene_list_fc <- function(input_file, species) {
  org  <- get_org(species)

  # always temp; clean afterwards
  tdir <- tempfile("kegg_"); dir.create(tdir)
  on.exit(unlink(tdir, recursive = TRUE, force = TRUE), add = TRUE)
  oldwd <- setwd(tdir); on.exit(setwd(oldwd), add = TRUE)

  deg <- read.table(input_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
  if (!all(c("Gene.ID","log2FC") %in% names(deg))) {
    stop("Input TXT must have columns: Gene.ID and log2FC")
  }
  deg <- deg |> distinct() |> na.omit()

  m  <- bitr(deg$Gene.ID, fromType="SYMBOL", toType="ENTREZID", OrgDb=org$orgdb)
  unmapped <- setdiff(deg$Gene.ID, m$SYMBOL)
  writeLines(unmapped, file.path(tdir,"unmapped_genes.txt"))
  mm <- merge(m, deg, by.x="SYMBOL", by.y="Gene.ID")
  if (!nrow(mm)) stop("No genes could be mapped to ENTREZ IDs.")

  ek <- tryCatch(enrichKEGG(gene = mm$ENTREZID, organism = org$kegg, pvalueCutoff = 0.05),
                 error = function(e) NULL)
  enr_df <- if (!is.null(ek)) as.data.frame(ek@result) else data.frame()
  full_csv <- file.path(tdir, "full_KEGG_enrichment_results.csv")
  write.csv(enr_df, full_csv, row.names = FALSE)

  barplot_path <- NULL
  top_csv      <- NULL

  if (nrow(enr_df) > 0) {
    top10 <- enr_df |> arrange(p.adjust) |> head(10)
    top_csv <- file.path(tdir, "top10_KEGG_enrichment.csv")
    write.csv(top10, top_csv, row.names = FALSE)

    top10$Description <- factor(top10$Description, levels = rev(top10$Description))
    p <- ggplot(top10, aes(x = Count, y = Description, fill = p.adjust)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "#d73027", high = "#4575b4", name = "p.adjust") +
      labs(title = "Top 10 KEGG Pathways", x = "Gene Count", y = NULL) +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            axis.text.y = element_text(size = 12))
    barplot_path <- file.path(tdir, "KEGG_Top10_Barplot.png")
    ggsave(barplot_path, p, width = 10, height = 6, dpi = 300, bg="white")

    # vectors for maps
    highlight <- setNames(rep(1, length(mm$ENTREZID)), mm$ENTREZID)
    gene_fc   <- setNames(as.numeric(mm$log2FC), mm$ENTREZID)
    gene_fc   <- gene_fc[is.finite(gene_fc)]

    for (i in seq_len(nrow(top10))) {
      pid  <- gsub(org$kegg, "", top10$ID[i])
      desc <- as.character(top10$Description[i])
      safe <- gsub("[^A-Za-z0-9_\\-]", "_", desc)

      # base highlight map
      try(pathview(gene.data  = highlight,
                   pathway.id = pid,
                   species    = org$kegg,
                   out.suffix = paste0(safe, "_BASE"),
                   kegg.native= TRUE,
                   same.layer = TRUE,
                   gene.idtype= "entrez",
                   limit      = list(gene = 1, cpd = 1),
                   low  = list(gene = NA),
                   mid  = list(gene = NA),
                   high = list(gene = NA),
                   kegg.dir   = tdir), silent = TRUE)

      # fold change colored map
      try(pathview(gene.data  = gene_fc,
                   pathway.id = pid,
                   species    = org$kegg,
                   out.suffix = paste0(safe, "_FC"),
                   kegg.native= TRUE,
                   same.layer = TRUE,
                   gene.idtype= "entrez",
                   kegg.dir   = tdir), silent = TRUE)
    }
  }

  # collect maps
  all_png <- list.files(tdir, pattern="\\.png$", full.names = TRUE)
  base_maps <- all_png[ is_base(all_png, org$kegg) & !is_fc(all_png) ]
  fc_maps   <- all_png[ is_fc(all_png) ]

  # zip
  zip_path <- file.path(tdir, "kegg_results.zip")
  files_to_zip <- c(full_csv,
                    if (length(top_csv)) top_csv,
                    file.path(tdir,"unmapped_genes.txt"),
                    if (length(barplot_path)) barplot_path,
                    base_maps, fc_maps)
  flat_zip(zip_path, files_to_zip, root = tdir)

  list(
    success             = TRUE,
    mode                = "genelist_fc",
    species             = species,
    barplot_base64      = if (length(barplot_path)) b64(barplot_path) else NULL,
    maps                = setNames(lapply(base_maps, b64), basename(base_maps)),
    foldchange_maps     = setNames(lapply(fc_maps,   b64), basename(fc_maps)),
    zip_base64          = b64(zip_path),
    zip_filename        = "kegg_genelist_fc_results.zip"
  )
}
