
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
    list(orgdb = org.Mm.eg.db, kegg = "mmu")
  } else if (species == "human") {
    list(orgdb = org.Hs.eg.db, kegg = "hsa")
  } else {
    stop("Unsupported species: must be 'mouse' or 'human'")
  }
}

# ===============================================================
# Function: Run KEGG from gene list only (highlight only)
# ===============================================================
run_kegg_from_gene_list <- function(input_file, species) {
  output_dir <- "Keggmap_output"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  org <- get_org(species)

  # read gene symbols
  gene_symbols <- read_lines(input_file) %>% unique() %>% na.omit()
  gene_df <- bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org$orgdb
  )

  # write unmapped
  unmapped <- setdiff(gene_symbols, gene_df$SYMBOL)
  writeLines(unmapped, file.path(output_dir, "unmapped_genes.txt"))

  # run enrichment
  kegg_enrich <- enrichKEGG(
    gene = gene_df$ENTREZID,
    organism = org$kegg,
    pvalueCutoff = 0.05
  )

  enrich_df <- as.data.frame(kegg_enrich@result)
  write.csv(enrich_df,
            file = file.path(output_dir, "full_KEGG_enrichment_results.csv"),
            row.names = FALSE)

  barplot_path <- NULL
  top10_csv <- NULL
  map_files <- character(0)
  map_files_base64 <- list()

  if (nrow(enrich_df) > 0) {
    # create top10
    top10 <- enrich_df %>% arrange(p.adjust) %>% head(10)
    top10_csv <- file.path(output_dir, "top10_KEGG_enrichment.csv")
    write.csv(top10, top10_csv, row.names = FALSE)

    # barplot
    top10$Description <- factor(top10$Description, levels = rev(top10$Description))
    barplot <- ggplot(top10, aes(x = Count, y = Description, fill = p.adjust)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "#d73027", high = "#4575b4", name = "p.adjust") +
      labs(title = "Top 10 KEGG Pathways", x = "Gene Count", y = NULL) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 12)
      )
    barplot_path <- file.path(output_dir, "KEGG_Top10_Barplot.png")
    ggsave(filename = barplot_path, plot = barplot, width = 10, height = 6)

    # prepare highlight-only gene data
    highlight_data <- setNames(rep(1, length(gene_df$ENTREZID)), gene_df$ENTREZID)

    # generate pathway maps
    for (i in seq_len(nrow(top10))) {
      pid <- top10$ID[i]
      tryCatch({
        pathview(
          gene.data = highlight_data,
          pathway.id = gsub(org$kegg, "", pid),
          species = org$kegg,
          # out.suffix = paste0("KEGG_", desc),
          kegg.native = TRUE,
          same.layer = TRUE,
          gene.idtype = "entrez",
          limit = list(gene = 1, cpd = 1),
          low = list(gene = NA),
          mid = list(gene = NA),
          high = list(gene = NA),
          kegg.dir = output_dir
        )
      }, error = function(e) {
        message("Map error: ", e$message)
      })
    }


# collect only clean KEGG maps (explicitly exclude foldchange in case)
all_pngs <- list.files(output_dir, pattern = "\\.png$", full.names = TRUE)
map_files <- all_pngs[
  grepl(paste0("^", org$kegg, "[0-9]+.*\\.png$"), basename(all_pngs)) &
  !grepl("foldchange", basename(all_pngs), ignore.case = TRUE)
]

# convert each map to base64
map_files_base64 <- list()
if (length(map_files) > 0) {
  for (f in map_files) {
    if (file.exists(f)) {
      img_bytes <- readBin(f, what = "raw", n = file.info(f)$size)
      map_files_base64[[basename(f)]] <- base64encode(img_bytes)
    }
  }
}

  }

  # barplot base64
  barplot_base64 <- NULL
  if (!is.null(barplot_path) && file.exists(barplot_path)) {
    img_bytes <- readBin(barplot_path, what = "raw", n = file.info(barplot_path)$size)
    barplot_base64 <- base64encode(img_bytes)
  }

  # zip everything
  zip_file <- file.path(output_dir, "kegg_results.zip")
  if (file.exists(zip_file)) file.remove(zip_file)
  zip::zip(zipfile = zip_file, files = list.files(output_dir, full.names = TRUE))

  return(list(
    barplot_base64 = barplot_base64,
    top10_csv = top10_csv,
    full_csv = file.path(output_dir, "full_KEGG_enrichment_results.csv"),
    unmapped = file.path(output_dir, "unmapped_genes.txt"),
    maps = map_files_base64,
    foldchange_maps = character(0),
    zip_file = zip_file,
    output_dir = output_dir
  ))
}

# ===============================================================
# Function 2: Run KEGG from gene list + fold change
# ===============================================================
run_kegg_from_gene_list_fc <- function(input_file, species) {
  output_dir <- file.path(getwd(), "Keggmap_output")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    }

    org <- get_org(species)

  # 📥 Read a tab‑delimited TXT with columns: Gene.ID and log2FC
    deg_data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # 🧹 Clean up the data
    if (!all(c("Gene.ID", "log2FC") %in% colnames(deg_data))) {
      stop("Input TXT must have columns: Gene.ID and log2FC")
    }
    deg_data <- deg_data %>% na.omit() %>% unique()
    # 🔎 Map SYMBOL → ENTREZ
    gene_df <- bitr(deg_data$Gene.ID,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org$orgdb)

    unmapped_genes <- setdiff(deg_data$Gene.ID, gene_df$SYMBOL)
    writeLines(unmapped_genes, file.path(output_dir, "unmapped_genes.txt"))

    # 🔀 Merge fold-change data with mapped IDs
    deg_merged <- merge(gene_df, deg_data, by.x = "SYMBOL", by.y = "Gene.ID")


    kegg_enrich <- enrichKEGG(gene = deg_merged$ENTREZID, organism = org$kegg, pvalueCutoff = 0.05)

    barplot_path <- NULL
    top10_csv <- NULL
    map_files <- character(0)
    fc_map_files <- character(0)

    if (!is.null(kegg_enrich@result) && nrow(kegg_enrich@result) > 0) {
      write.csv(as.data.frame(kegg_enrich@result),
                file = file.path(output_dir, "full_KEGG_enrichment_results.csv"),
                row.names = FALSE)

      kegg_top10 <- head(kegg_enrich@result[order(kegg_enrich@result$p.adjust), ], 10)
      top10_csv <- file.path(output_dir, "top10_KEGG_enrichment.csv")
      write.csv(kegg_top10, top10_csv, row.names = FALSE)

      kegg_top10$Description <- factor(kegg_top10$Description, levels = rev(kegg_top10$Description))
      barplot_path <- file.path(output_dir, "KEGG_Top10_Barplot.png")
      cat("barplot_path: ", barplot_path, "\n")
      png(filename = barplot_path, width = 1200, height = 800, res = 150)

      print(
        ggplot(kegg_top10, aes(x = Description, y = -log10(p.adjust))) +
          geom_bar(stat = "identity", fill = "#1F77B4") +
          coord_flip() +
          labs(title = "Top 10 Enriched KEGG Pathways",
              x = "Pathway",
              y = expression(-log[10]~"(adjusted p-value)")) +
          theme_minimal(base_size = 14)
      )
      dev.off()

      # After saving the plot, read the image as raw bytes and return as base64 string
      img_bytes <- readBin(barplot_path, what = "raw", n = file.info(barplot_path)$size)
      barplot_base64 <- base64enc::base64encode(img_bytes)


      gene_fc <- setNames(as.numeric(deg_merged$log2FC), deg_merged$ENTREZID)
      gene_fc <- gene_fc[!is.na(gene_fc) & is.finite(gene_fc)]
    
      # Save current working directory
      old_wd <- getwd()

      # Change working directory to your output folder
      setwd(output_dir)
      for (i in seq_len(nrow(kegg_top10))) {
        pid <- kegg_top10$ID[i]
        desc <- kegg_top10$Description[i]

        # Sanitize and append "_FC"
        safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", desc)
        message("🔧 Generating map for: ", pid, " with suffix: ", safe_name)

        tryCatch({
          pathview(
            gene.data = gene_fc,
            pathway.id = pid,
            species = org$kegg,
            out.suffix = safe_name,
            gene.idtype = "entrez",
            limit = list(gene = 3),
            kegg.dir = output_dir,
          )
        }, error = function(e) {
          message("Map error: ", e$message, "\n")
        })
      }



      setwd(old_wd)

      # collect maps
      all_maps <- list.files(output_dir, pattern = "\\.png$", full.names = TRUE)
      map_files <- all_maps[grepl(paste0("^", org$kegg, "[0-9]+.*\\.png$"), basename(all_maps))]
      
      cat("==============================================\n")
      cat("all_maps: ", all_maps, "\n")
      cat("map_files: ", map_files, "\n")
      cat("==============================================\n")
      fc_map_files <- all_maps[grepl("foldchange", all_maps, ignore.case = TRUE)]
      # Convert map_files (list of PNGs) to base64
      map_files_base64 <- list()
      if (length(map_files) > 0) {
        for (f in map_files) {
          if (file.exists(f)) {
            img_bytes <- readBin(f, what = "raw", n = file.info(f)$size)
            img_base64 <- base64enc::base64encode(img_bytes)
            map_files_base64[[f]] <- img_base64
          } else {
            map_files_base64[[f]] <- NULL
            message("File does not exist: ", f)
          }
        }
      }

      # Convert fc_map_files (list of PNGs) to base64
      fc_map_files_base64 <- list()
      if (length(fc_map_files) > 0) {
        for (f in fc_map_files) {
          if (file.exists(f)) {
            img_bytes <- readBin(f, what = "raw", n = file.info(f)$size)
            img_base64 <- base64enc::base64encode(img_bytes)
            fc_map_files_base64[[f]] <- img_base64
          } else {
            fc_map_files_base64[[f]] <- NULL
            message("File does not exist: ", f)
          }
        }
      }
    }

    # zip
    zip_file <- file.path(output_dir, "kegg_results.zip")
    if (file.exists(zip_file)) file.remove(zip_file)
    zip::zip(zipfile = zip_file, files = list.files(output_dir, full.names = TRUE))

    return(list(
      barplot_base64 = barplot_base64,
      top10_csv = top10_csv,
      full_csv = file.path(output_dir, "full_KEGG_enrichment_results.csv"),
      unmapped = file.path(output_dir, "unmapped_genes.txt"),
      maps = map_files_base64,
      foldchange_maps = fc_map_files_base64,
      zip_file = zip_file,
      output_dir = output_dir
    ))
}
