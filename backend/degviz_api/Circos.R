# Circos.R
suppressPackageStartupMessages({
  library(plumber)
  library(jsonlite)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(circlize)
  library(base64enc)
  library(tools)
})

# ---------- utils ----------

`%||%` <- function(a, b) if (!is.null(a)) a else b

trySuppress <- function(expr) suppressWarnings(try(expr, silent = TRUE))

normalize_genes <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\.\\d+$", "", x)   # drop Ensembl version suffix
  toupper(x)                     # unify case
  x
}

read_any <- function(path) {
  stopifnot(file.exists(path))
  ext <- tolower(file_ext(path))
  if (ext %in% c("csv")) {
    df <- trySuppress(read_csv(path, show_col_types = FALSE))
  } else if (ext %in% c("tsv", "txt")) {
    df <- trySuppress(read_tsv(path, show_col_types = FALSE))
    if (inherits(df, "try-error")) df <- trySuppress(read_csv(path, show_col_types = FALSE))
    if (inherits(df, "try-error")) df <- read_table2(path, show_col_types = FALSE)
  } else {
    df <- trySuppress(read_tsv(path, show_col_types = FALSE))
    if (inherits(df, "try-error")) df <- trySuppress(read_csv(path, show_col_types = FALSE))
    if (inherits(df, "try-error")) df <- read_table2(path, show_col_types = FALSE)
  }
  df
}

guess_gene_col <- function(df, user_gene_col = NULL) {
  # exact match if user provided and exists
  if (!is.null(user_gene_col) && user_gene_col %in% names(df)) return(user_gene_col)
  # try case-insensitive match
  if (!is.null(user_gene_col)) {
    nm <- names(df)
    hit <- nm[tolower(nm) == tolower(user_gene_col)]
    if (length(hit) >= 1) return(hit[[1]])
  }
  # common candidates (case-insensitive)
  candidates <- c("Gene","Gene.ID","gene","symbol","Symbol","GeneSymbol","GeneName","ENSEMBL","Ensembl","id","ID")
  nm <- names(df)
  hit <- nm[tolower(nm) %in% tolower(candidates)]
  if (length(hit) >= 1) return(hit[[1]])
  stop("Could not find a gene column in the uploaded file(s).")
}

build_binary_matrix <- function(named_gene_lists) {
  all_genes <- sort(unique(unlist(named_gene_lists)))
  if (!length(all_genes)) return(data.frame(Gene = character(0)))
  mat <- matrix(0L, nrow = length(all_genes), ncol = length(named_gene_lists),
                dimnames = list(all_genes, names(named_gene_lists)))
  for (lab in names(named_gene_lists)) {
    present <- intersect(all_genes, named_gene_lists[[lab]])
    if (length(present)) mat[present, lab] <- 1L
  }
  as.data.frame(mat) |> rownames_to_column("Gene")
}

build_links_from_matrix <- function(binary_df, min_shared = 1L) {
  if (!nrow(binary_df)) return(data.frame(from=character(), to=character(), value=integer()))
  bm <- binary_df
  rownames(bm) <- bm$Gene; bm$Gene <- NULL
  m <- as.matrix(bm); labs <- colnames(m)
  if (ncol(m) < 2) return(data.frame(from=character(), to=character(), value=integer()))
  out <- list(); k <- 1L
  for (i in 1:(ncol(m)-1)) for (j in (i+1):ncol(m)) {
    v <- sum(m[, i] == 1 & m[, j] == 1)
    if (v >= min_shared) {
      out[[k]] <- data.frame(from=labs[i], to=labs[j], value=v, stringsAsFactors = FALSE); k <- k + 1L
    }
  }
  if (!length(out)) return(data.frame(from=character(), to=character(), value=integer()))
  bind_rows(out)
}

pairwise_overlap <- function(gene_lists) {
  labs <- names(gene_lists)
  if (length(labs) < 2) return(matrix(0, 0, 0))
  m <- matrix(0L, nrow = length(labs), ncol = length(labs), dimnames = list(labs, labs))
  for (i in seq_along(labs)) for (j in seq_along(labs)) {
    if (i <= j) m[i, j] <- length(intersect(gene_lists[[i]], gene_lists[[j]])) else m[i, j] <- m[j, i]
  }
  m
}

render_chord_base64 <- function(links_df) {
  if (nrow(links_df) == 0) stop("No links to plot.")
  all_sectors <- unique(c(links_df$from, links_df$to))

  pdf_path <- tempfile(fileext = ".pdf")
  png_path <- tempfile(fileext = ".png")

  # PDF
  grDevices::pdf(pdf_path, width = 8, height = 8, onefile = TRUE)
  circos.clear()
  lp <- links_df
  lp$from <- factor(lp$from, levels = all_sectors)
  lp$to   <- factor(lp$to,   levels = all_sectors)
  chordDiagram(lp[, c("from","to","value")], annotationTrack="grid", preAllocateTracks=1)
  circos.trackPlotRegion(track.index=1, panel.fun=function(x,y){
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing="clockwise", niceFacing=TRUE, adj=c(0,0.5), cex=0.6)
  }, bg.border=NA)
  grDevices::dev.off()

  # PNG
  grDevices::png(png_path, width = 1600, height = 1600, res = 200)
  circos.clear()
  lp <- links_df
  lp$from <- factor(lp$from, levels = all_sectors)
  lp$to   <- factor(lp$to,   levels = all_sectors)
  chordDiagram(lp[, c("from","to","value")], annotationTrack="grid", preAllocateTracks=1)
  circos.trackPlotRegion(track.index=1, panel.fun=function(x,y){
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing="clockwise", niceFacing=TRUE, adj=c(0,0.5), cex=0.6)
  }, bg.border=NA)
  grDevices::dev.off()

  on.exit({ unlink(pdf_path, force = TRUE); unlink(png_path, force = TRUE) }, add = TRUE)
  pdf_b <- readBin(pdf_path, "raw", n = file.info(pdf_path)$size)
  png_b <- readBin(png_path, "raw", n = file.info(png_path)$size)

  list(pdf_base64 = base64encode(pdf_b), png_base64 = base64encode(png_b))
}

