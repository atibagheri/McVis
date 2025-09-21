# plumber.R
# Enable CORS

library(plumber)
library(jsonlite)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(VennDiagram)
library(UpSetR)
library(futile.logger)
library(pathview)
library(ggplot2)
library(readr)
library(dplyr)
library(openxlsx)
library(zip)
library(rentrez)
library(httr)
library(tibble)
library(stringr)
library(circlize)
library(tools)
library(base64enc)

source("Transpca.R")
source("Venn.R")
source("TextMining.R")
source("Gomap.R")
source("Keggmap.R")
source("Circos.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b


#* @filter cors
function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  res$setHeader("Access-Control-Allow-Methods", "GET, POST, OPTIONS, PUT, DELETE")
  res$setHeader("Access-Control-Allow-Headers", "Content-Type, Authorization")

  if (req$REQUEST_METHOD == "OPTIONS") {
    res$status <- 200
    return(list())
  } else {
    plumber::forward()
  }
}



########################################################
#Mouse Macrophage Polarization
###########################################
# ---------- ROUTES ----------

#* @get /healthz
#* @serializer unboxedJSON
function() list(status="ok", build="transpca-all-in-one")

# 5-row preview (array of row objects)
#* @get /mouse-tpm-preview
#* @serializer unboxedJSON
function(limit = 5) {
  can <- canonicalize_table_path(DEFAULT_MOUSE_TPM)
  if (!can$ok) return(list(status="error", message=sprintf("Mouse TPM file not found at: %s", DEFAULT_MOUSE_TPM)))
  lim <- suppressWarnings(as.integer(limit[1])); if (is.na(lim) || lim < 1) lim <- 5
  prev <- tryCatch(read_table_preview(can$path, lim), error = function(e) NULL)
  if (is.null(prev)) return(list(status="error", message="Failed to read mouse TPM file"))
  rows <- lapply(seq_len(nrow(prev)), function(i) as.list(prev[i, , drop = FALSE]))
  list(status="ok", columns=names(prev), rows=rows)
}

# Full download (raw bytes, text/plain; suggest .txt filename)
#* @get /mouse-tpm-download
#* @serializer contentType list(type="text/plain; charset=UTF-8")
function(res) {
  path <- DEFAULT_MOUSE_TPM
  if (!file.exists(path)) { res$status <- 404; return("Not found") }
  bn <- basename(path)
  if (!grepl("\\.txt$", bn, ignore.case = TRUE)) {
    bn <- paste0(tools::file_path_sans_ext(bn), ".txt")
  }
  res$setHeader("Content-Disposition", sprintf('attachment; filename="%s"', bn))
  readBin(path, what = "raw", n = file.info(path)$size)
}

# JSON upload (base64 fields). Accepts:
# {
#   "file": {"name": "...", "b64": "data:...;base64,...."},
#   "meta": {"name": "...", "b64": "data:...;base64,...."},
#   "mouse": {"name": "...", "b64": "..."}  (optional, overrides default)
#   "size_col": "...", "color_col": "...", "shape_col": "...", "label_col": "..."
# }
#* @post /transpca_json
#* @serializer unboxedJSON
function(req, res) {
  parsed <- tryCatch(jsonlite::fromJSON(req$postBody, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(parsed)) { res$status <- 400; return(list(status="error", message="Invalid JSON body")) }

  strip_data_uri <- function(s) sub("^data:.*;base64,", "", s %||% "")
  write_b64_to_temp <- function(name, b64) {
    if (is.null(b64) || !nzchar(b64)) return(NULL)
    b64 <- strip_data_uri(b64)
    ext <- tolower(tools::file_ext(name %||% "txt"))
    suffix <- paste0(".", ifelse(nzchar(ext), ext, "txt"))
    tmp <- tempfile(fileext = suffix)
    raw <- tryCatch(base64enc::base64decode(b64), error = function(e) NULL)
    if (is.null(raw)) return(NULL)
    writeBin(raw, tmp); tmp
  }

  user_path  <- write_b64_to_temp(parsed$file$name  %||% parsed$file_name  %||% "user.txt",
                                  parsed$file$b64   %||% parsed$file_b64)
  meta_path  <- write_b64_to_temp(parsed$meta$name  %||% parsed$meta_name,
                                  parsed$meta$b64   %||% parsed$meta_b64)
  mouse_path <- write_b64_to_temp(parsed$mouse$name %||% parsed$mouse_name,
                                  parsed$mouse$b64  %||% parsed$mouse_b64) %||% DEFAULT_MOUSE_TPM

  if (is.null(user_path)) { res$status <- 400; return(list(status="error", message="Missing or invalid 'file'")) }

  tryCatch({
    run_transpca_core(
      user_path,
      mouse_tpm = mouse_path,
      meta_path  = meta_path,
      size_col   = parsed$size_col %||% NULL,
      color_col  = parsed$color_col %||% NULL,
      shape_col  = parsed$shape_col %||% NULL,
      label_col  = parsed$label_col %||% NULL
    )
  }, error = function(e) { res$status <- 400; list(status="error", message=as.character(e)) })
}

# Absolute-path POST (for Flask proxy). Fields:
#   file_path=/ABS/PATH/matrix.txt
#   meta_path=/ABS/PATH/meta.txt   (optional)
#   size_col=Day (optional)
#   color_col=WoundZone (optional)
#   shape_col=...  label_col=...
#* @post /transpca
#* @serializer unboxedJSON
function(file_path = "", meta_path = "", size_col = "", color_col = "", shape_col = "", label_col = "") {
  if (!nzchar(file_path)) {
    return(list(
      status  = "error",
      message = "Provide absolute paths (no uploads). Example:\n  curl -sS http://127.0.0.1:8000/transpca \\\n    -d file_path=/ABS/PATH/Canine.txt \\\n    -d meta_path=/ABS/PATH/vml_metadata.txt \\\n    -d label_col=SampleName"
    ))
  }
  tryCatch({
    run_transpca_core(
      file_path,
      mouse_tpm = DEFAULT_MOUSE_TPM,
      meta_path = if (nzchar(meta_path)) meta_path else NULL,
      size_col  = if (nzchar(size_col))  size_col  else NULL,
      color_col = if (nzchar(color_col)) color_col else NULL,
      shape_col = if (nzchar(shape_col)) shape_col else NULL,
      label_col = if (nzchar(label_col)) label_col else NULL
    )
  }, error = function(e) list(status="error", message=as.character(e)))
}

########################################################
#Venn/UpSet
########################################################


#* @post /venn-upset
function(req, res) {
  cat("==== /venn-upset CALLED ====\n")
  parsed <- tryCatch({ jsonlite::fromJSON(req$postBody) }, error = function(e) NULL)
  if (is.null(parsed) || is.null(parsed$paths) || length(parsed$paths) < 2) {
    res$status <- 400
    return(list(error = "At least two file paths required"))
  }

  output_path_png <- tempfile(fileext = ".png")
  output_path_pdf <- tempfile(fileext = ".pdf")
  # open the first file to inspect
  result <- generate_venn_or_upset(parsed$paths, output_path_png, output_path_pdf)

  # After saving both PNG and PDF
  pdf_bytes <- readBin(output_path_pdf, what = "raw", n = file.info(output_path_pdf)$size)
  pdf_base64 <- base64enc::base64encode(pdf_bytes)
  png_bytes <- readBin(output_path_png, what = "raw", n = file.info(output_path_png)$size)
  png_base64 <- base64enc::base64encode(png_bytes)

return(list(
  success = TRUE,
  output_png = result$output_png,
  output_pdf = result$output_pdf,
  png_base64 = result$png_base64,
  pdf_base64 = pdf_base64
  ))
}


########################################################
#GOMAP 
########################################################

#* @post /gomap
#* @param mode The mode: "genelist" or "genelist_fc"
#* @param species The species: "mouse" or "human"
#* @param file_path Path to the input gene file
function(req, res) {
  cat("==== /gomap CALLED ====\n")
  parsed <- tryCatch(jsonlite::fromJSON(req$postBody), error = function(e) NULL)
  if (is.null(parsed) || is.null(parsed$file_path) || is.null(parsed$mode) || is.null(parsed$species)) {
    res$status <- 400
    return(list(success = FALSE, error = "Missing required parameters: file_path, mode, species"))
  }

  mode <- parsed$mode; species <- parsed$species; file <- parsed$file_path
  cat("📄 file_path:", file, "\n"); cat("🔧 mode:", mode, "\n"); cat("🧬 species:", species, "\n")

  result <- tryCatch({
    if (identical(mode, "genelist"))      run_go_from_gene_list(file, species)
    else if (identical(mode, "genelist_fc")) run_go_from_gene_list_fc(file, species)
    else stop("Invalid mode. Use 'genelist' or 'genelist_fc'.")
  }, error = function(e) { res$status <- 500; return(list(success = FALSE, error = e$message)) })

  if (isFALSE(result$success %||% TRUE)) { res$status <- 500; return(result) }

  cat("✅ zip_base64 length:", nchar(result$zip_base64 %||% ""), "\n")

  # Whitelist ONLY safe fields (no absolute paths)
  list(
    success            = TRUE,
    mode               = mode,
    species            = species,
    n_input_genes      = result$n_input_genes,
    n_mapped           = result$n_mapped,
    n_terms            = result$n_terms,
    barplot_base64     = result$barplot_base64,
    cnetplot_base64    = result$cnetplot_base64,
    zip_base64         = result$zip_base64,        # <-- one ZIP with CSV+plots+unmapped
    zip_filename       = result$zip_filename,
    csv_base64         = result$csv_base64,        # (optional)
    csv_filename       = result$csv_filename,
    unmapped_base64    = result$unmapped_base64,   # (optional)
    unmapped_filename  = result$unmapped_filename
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x


########################################################
#KEGGMAP 
########################################################
# plumber_kegg.R

#* @apiTitle KEGG Map API
#* @serializer unboxedJSON

#* Run KEGG
#* @post /keggmap
function(req, res) {
  parsed <- tryCatch(jsonlite::fromJSON(req$postBody), error = function(e) NULL)
  if (is.null(parsed) || is.null(parsed$file_path) || is.null(parsed$mode) || is.null(parsed$species)) {
    res$status <- 400
    return(list(success = FALSE, error = "Missing required parameters: file_path, mode, species"))
  }

  mode    <- parsed$mode
  species <- parsed$species
  file    <- parsed$file_path

  out <- tryCatch({
    if (identical(mode, "genelist")) {
      run_kegg_from_gene_list(file, species)
    } else if (identical(mode, "genelist_fc")) {
      run_kegg_from_gene_list_fc(file, species)
    } else {
      stop("Invalid mode. Use 'genelist' or 'genelist_fc'.")
    }
  }, error = function(e) {
    res$status <- 500
    return(list(success = FALSE, error = e$message))
  })

  # Unbox small scalars explicitly (if your plumber doesn't auto_unbox)
  unbox <- jsonlite::unbox
  list(
    success            = TRUE,
    mode               = unbox(mode),
    species            = unbox(species),
    n_input_genes      = unbox(out$n_input_genes),
    n_mapped           = unbox(out$n_mapped),
    n_terms            = unbox(out$n_terms),
    barplot_base64     = out$barplot_base64,
    csv_base64         = out$csv_base64,
    csv_filename       = out$csv_filename,
    top10_csv_base64   = out$top10_csv_base64,
    top10_csv_filename = out$top10_csv_filename,
    unmapped_base64    = out$unmapped_base64,
    unmapped_filename  = out$unmapped_filename,
    maps               = out$maps,
    foldchange_maps    = out$foldchange_maps,
    zip_base64         = out$zip_base64,
    zip_filename       = out$zip_filename
  )
}

########################################################
#Text Mining
########################################################

# ---------- Health check ----------
#* Health check
#* @get /health
function() {
  list(ok = TRUE, time = as.character(Sys.time()))
}

# ---------- Main endpoint ----------
#* PubMed text mining (genes file + optional term)
#*
#* @serializer json list(na="null", auto_unbox=TRUE, pretty=FALSE)
#* @post /api/textmining
function(req, res) {
  # Multipart form fields:
  # - genes: file (required)
  # - mode: "gene_only" | "gene_term" (optional, default gene_only)
  # - term: string (optional when mode=gene_term)

  # Validate upload
  if (is.null(req$files) || is.null(req$files$genes) || is.na(req$files$genes$datapath)) {
    res$status <- 400
    return(list(success = FALSE, error = "Missing 'genes' file (multipart/form-data)."))
  }

  genes_path <- req$files$genes$datapath
  mode <- req$body$mode %||% "gene_only"
  term <- req$body$term %||% ""

  # Sanity checks
  if (!file.exists(genes_path)) {
    res$status <- 400
    return(list(success = FALSE, error = "Uploaded file not found on server."))
  }
  if (!(mode %in% c("gene_only", "gene_term"))) {
    res$status <- 400
    return(list(success = FALSE, error = "Invalid 'mode'. Use 'gene_only' or 'gene_term'."))
  }

  # Run analysis (wrap in tryCatch to return JSON errors)
  out <- tryCatch({
    run_pubmed_analysis(
      file_path = genes_path,
      mode      = mode,
      term      = term,
      # You can tune these defaults if needed:
      top_n_per_gene          = 3,
      sleep_sec               = 0.34,
      wc_seed                 = 42,
      wc_max_words            = 200,
      wc_min_freq             = 1,
      wc_scale                = c(4, 0.7),
      summaries_for_zero_hits = FALSE
    )
  }, error = function(e) {
    res$status <- 500
    list(success = FALSE, error = paste("Server error:", conditionMessage(e)))
  })

  # Ensure response is always JSON serializable
  out
}

# helper for null-coalescing (keeps this file self-contained)
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b


########################################################
#Circos
#######################################################
#* Generate circos chord plot from multiple gene lists (no persistent outputs)
#* @post /circos
#* @serializer json
function(req, res) {
  cat("==== /Circos CALLED ====\n")
  parsed <- tryCatch(fromJSON(req$postBody, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(parsed) || is.null(parsed$paths) || length(parsed$paths) < 2) {
    res$status <- 400; return(list(error = "At least two file paths required"))
  }

  paths         <- parsed$paths
  gene_col_user <- parsed$gene_col_user %||% parsed$gene_col %||% NULL
  min_shared    <- if (!is.null(parsed$min_shared)) as.integer(parsed$min_shared) else 1L

  labels <- parsed$labels
  if (is.null(labels) || length(labels) != length(paths)) {
    labels <- vapply(paths, function(x) sub("\\.[^.]+$", "", basename(x)), character(1))
  }
  labels <- as.character(labels)

  # build gene lists
  gene_lists <- setNames(vector("list", length(paths)), labels)
  used_cols <- character(length(paths))
  for (i in seq_along(paths)) {
    df <- read_any(paths[i])
    gc <- guess_gene_col(df, gene_col_user)
    used_cols[i] <- gc
    genes <- normalize_genes(df[[gc]]) |> unique() |> stats::na.omit()
    gene_lists[[ labels[i] ]] <- genes
  }

  # matrix + links
  binary_df <- build_binary_matrix(gene_lists)
  links_df  <- build_links_from_matrix(binary_df, as.integer(min_shared))

  if (!nrow(links_df)) {
    res$status <- 422
    return(list(
      error = "No links to plot",
      hint  = "Check gene column and make sure lists use the same ID system (SYMBOL vs ENSEMBL).",
      gene_column_used = as.list(stats::setNames(as.list(used_cols), labels)),
      sizes = as.list(vapply(gene_lists, length, 1L)),
      pairwise_overlap = pairwise_overlap(gene_lists)
    ))
  }

  out <- render_chord_base64(links_df)
  list(success = TRUE, png_base64 = out$png_base64, pdf_base64 = out$pdf_base64)
}