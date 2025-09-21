# TextMining.R — complete, polished, and matched to your React frontend
# - Shows TOP N genes (default 20) in plots for readability
# - CSVs include ALL genes
# - Summaries include: Gene, PMID, Title, Journal, Year, Score, URL
# - Returns JSON-ready list with keys React expects

suppressPackageStartupMessages({
  library(rentrez)
  library(ggplot2)
  library(base64enc)
  library(wordcloud)
  library(RColorBrewer)
  library(dplyr)
  library(stringr)
  library(purrr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# (optional) NCBI etiquette:
# Sys.setenv(ENTREZ_EMAIL = "you@example.org")
# Sys.setenv(ENTREZ_KEY   = "<your_api_key>")
if (nzchar(Sys.getenv("ENTREZ_KEY", ""))) entrez_key(Sys.getenv("ENTREZ_KEY"))

# ---------- Helpers ----------
safe_entrez_search <- function(term, retmax = 0, tries = 3, sleep_sec = 0.34) {
  for (k in seq_len(tries)) {
    out <- tryCatch(entrez_search(db = "pubmed", term = term, retmax = retmax),
                    error = function(e) NULL)
    if (!is.null(out)) return(out)
    Sys.sleep(sleep_sec * k) # simple backoff
  }
  NULL
}

# Summaries with requested fields (+ URL)
get_pubmed_info <- function(gene, extra_term = "", max_results = 5) {
  tryCatch({
    query <- if (nzchar(extra_term)) paste(gene, extra_term, sep = " AND ") else gene
    res <- entrez_search(db = "pubmed", term = query, retmax = max_results)
    ids <- res$ids
    if (length(ids) == 0) return(NULL)

    summaries <- entrez_summary(db = "pubmed", id = ids)
    rows <- lapply(seq_along(summaries), function(i) {
      sm <- summaries[[i]]
      data.frame(
        Gene    = gene,
        PMID    = sm$uid %||% NA,
        Title   = sm$title %||% NA,
        Year    = sm$pubdate %||% NA,      # not parsed (as requested)
        Score   = max_results - i + 1,
        stringsAsFactors = FALSE
      )
    })
    out <- bind_rows(rows)
    out$URL <- paste0("https://pubmed.ncbi.nlm.nih.gov/", out$PMID, "/")
    out
  }, error = function(e) {
    message("Error for gene ", gene, ": ", e$message)
    NULL
  })
}

# ---------- Main ----------
run_pubmed_analysis <- function(
  file_path,
  mode                   = "gene_only",   # "gene_only" | "gene_term"
  term                   = "",
  top_n_per_gene         = 3,
  top_n_genes            = 20,            # <-- NEW: plots show top N by hits
  sleep_sec              = 0.34,
  wc_seed                = 42,
  wc_max_words           = 200,
  wc_min_freq            = 1,
  wc_scale               = c(4, 0.7),
  inches_per_gene        = 0.33,          # vertical space per gene (inches)
  min_plot_height        = 6,
  max_plot_height        = 40
) {
  # 0) Read + clean genes
  genes <- readLines(file_path, warn = FALSE)
  genes <- str_trim(genes)
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0) stop("Gene list is empty")

  # Normalize to upper + dedupe, preserve order
  g_upper <- toupper(genes)
  genes_unique <- g_upper[!duplicated(g_upper)]

  build_query <- function(g) {
    if (identical(mode, "gene_term") && nzchar(term)) paste0(g, " AND (", term, ")") else g
  }

  # 1) PubMed hit counts for all unique genes
  hit_counts_unique <- integer(length(genes_unique))
  for (i in seq_along(genes_unique)) {
    q <- build_query(genes_unique[i])
    res <- safe_entrez_search(q, retmax = 0, sleep_sec = sleep_sec)
    hit_counts_unique[i] <- if (is.null(res)) NA_integer_ else as.integer(res$count)
    Sys.sleep(sleep_sec)
  }
  hit_df <- tibble(
    Gene = genes_unique,
    PubMed_Hits = as.integer(replace_na(hit_counts_unique, 0L))
  )

  # 2) Subset for plotting/word cloud (top N genes)
  hit_df_plot <- hit_df %>%
    arrange(desc(PubMed_Hits)) %>%
    slice_head(n = top_n_genes)

  # 3) Dynamic plot height so labels never overlap
  n_genes <- nrow(hit_df_plot)
  plot_height <- max(min_plot_height, min(max_plot_height, n_genes * inches_per_gene))

  # 4) Bar plot (horizontal)
  bar_plot <- ggplot(hit_df_plot, aes(x = reorder(Gene, PubMed_Hits), y = PubMed_Hits)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.85) +
    coord_flip() +
    labs(
      title = if (mode == "gene_term" && nzchar(term))
        paste("Top", top_n_genes, "PubMed Hits: Gene AND", term)
      else
        paste("Top", top_n_genes, "Genes by PubMed Hits"),
      x = "Gene", y = "Number of Articles"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 8, hjust = 1),
      plot.margin = margin(t = 10, r = 20, b = 10, l = 90),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

  bar_png <- tempfile(fileext = ".png")
  bar_pdf <- tempfile(fileext = ".pdf")
  ggsave(bar_png, plot = bar_plot, width = 7, height = plot_height, dpi = 300)
  ggsave(bar_pdf, plot = bar_plot, width = 7, height = plot_height)

  # 5) Word cloud from the same top-N subset
  wc_dat <- hit_df_plot %>%
    filter(PubMed_Hits >= wc_min_freq) %>%
    arrange(desc(PubMed_Hits)) %>%
    slice_head(n = wc_max_words)

  set.seed(wc_seed)
  pal <- brewer.pal(8, "Dark2")

  wc_png <- tempfile(fileext = ".png")
  png(wc_png, width = 1000, height = 600, res = 150)
  par(mar = c(0, 0, 0, 0))
  wordcloud(
    words        = wc_dat$Gene,
    freq         = wc_dat$PubMed_Hits,
    min.freq     = wc_min_freq,
    random.order = FALSE,
    rot.per      = 0.15,
    colors       = pal,
    scale        = wc_scale
  )
  dev.off()

  wc_pdf <- tempfile(fileext = ".pdf")
  pdf(wc_pdf, width = 11, height = 8.5, useDingbats = FALSE)
  par(mar = c(0, 0, 0, 0))
  set.seed(wc_seed)
  wordcloud(
    words        = wc_dat$Gene,
    freq         = wc_dat$PubMed_Hits,
    min.freq     = wc_min_freq,
    random.order = FALSE,
    rot.per      = 0.15,
    colors       = pal,
    scale        = wc_scale
  )
  dev.off()

  # 6) Per‑gene summaries (ALL genes; requested fields)
  extra_term <- if (identical(mode, "gene_term") && nzchar(term)) term else ""
  summaries_list <- lapply(genes_unique, function(g) {
    Sys.sleep(sleep_sec)
    get_pubmed_info(g, extra_term = extra_term, max_results = top_n_per_gene)
  })
  summaries_df <- bind_rows(summaries_list)
  if (is.null(summaries_df)) {
    summaries_df <- tibble(
      Gene=character(), PMID=character(), Title=character(),
       Year=character(), Score=integer()
    )
  }

  # 7) CSVs (ALL genes)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  hits_csv      <- file.path(tempdir(), sprintf("pubmed_hits_%s.csv", ts))
  summaries_csv <- file.path(tempdir(), sprintf("pubmed_summaries_%s.csv", ts))
  write.csv(hit_df, hits_csv, row.names = FALSE)
  write.csv(summaries_df, summaries_csv, row.names = FALSE)

  # 8) Payload for React
  list(
    success = TRUE,
    params = list(
      mode = mode, term = term,
      top_n_per_gene = top_n_per_gene,
      top_n_genes = top_n_genes
    ),

    barplot_png_base64   = base64enc::base64encode(bar_png),
    barplot_pdf_base64   = base64enc::base64encode(bar_pdf),
    wordcloud_png_base64 = base64enc::base64encode(wc_png),
    wordcloud_pdf_base64 = base64enc::base64encode(wc_pdf),

    hitcount_csv_base64  = base64enc::base64encode(hits_csv),
    summary_csv_base64   = base64enc::base64encode(summaries_csv),

    table     = hit_df,        # ALL genes
    summaries = summaries_df,  # ALL genes

    downloads = list(          # URLs not used; React falls back to base64
      barplot_png   = NULL,
      barplot_pdf   = NULL,
      wordcloud_png = NULL,
      wordcloud_pdf = NULL,
      hitcount_csv  = NULL,
      summary_csv   = NULL
    )
  )
}
