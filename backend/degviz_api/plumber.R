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

# source("functions.R")
source("Venn.R")
source("Gomap.R")
source("Keggmap.R")


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
#Venn/UpSet
########################################################
source("Venn.R")

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
#ENRICHMENT
########################################################

#* @param file:file Gene list file (with column Gene.ID)
#* @param species The species: "mouse" or "human"
#* @post /enrichment
function(req, res) {
  cat("==== /enrichment CALLED ====\n")
  parsed <- tryCatch({ fromJSON(req$postBody) }, error = function(e) NULL)
  if (is.null(parsed) || is.null(parsed$file_path)) {
    res$status <- 400
    return(list(success = FALSE, error = "File path is missing"))
  }
  file_path <- parsed$file_path
  species <- parsed$species
  tryCatch({
    result <- run_enrichment_analysis(file_path, species)
    return(list(success = TRUE, 
    output_path = result$output_dir, 
    go_pdf = result$go_pdf, 
    kegg_pdf = result$kegg_pdf, 
    go_txt = result$go_txt, 
    kegg_txt = result$kegg_txt,
    go_base64 = result$go_base64,
    kegg_base64 = result$kegg_base64))
  }, error = function(e) {
    res$status <- 500
    return(list(success = FALSE, error = e$message))
  })
}



########################################################
#GOMAP 
########################################################
#* @param mode The mode: "genelist" or "genelist_fc"
#* @param species The species: "mouse" or "human"
#* @param file_path Path to the input gene file
#* @param output_dir Path to the output folder
#* @post /gomap
function(req, res) {
  cat("==== /gomap CALLED ====\n")
  
  # Parse JSON input
  parsed <- tryCatch({ jsonlite::fromJSON(req$postBody) }, error = function(e) NULL)
  
  # Validate required parameters
  if (is.null(parsed) || is.null(parsed$file_path) || is.null(parsed$mode) || is.null(parsed$species)) {
    res$status <- 400
    return(list(success = FALSE, error = "Missing required parameters: file_path, mode, species"))
  }
  
  mode <- parsed$mode
  species <- parsed$species
  file_path <- parsed$file_path
  
  cat("📄 file_path: ", file_path, "\n")
  cat("🔧 mode: ", mode, "\n")
  cat("🧬 species: ", species, "\n")
  
  # Run GO analysis
  tryCatch({
    if (mode == "genelist") {
      cat("==== Running in genelist mode ====\n")
      result <- run_go_from_gene_list(file_path, species)
    } else if (mode == "genelist_fc") {
      cat("==== Running in genelist_fc mode ====\n")
      result <- run_go_from_gene_list_fc(file_path, species)
    } else {
      stop("Invalid mode. Use 'genelist' or 'genelist_fc'.")
    }
    
    # ✅ Return results
    return(list(
      success = TRUE,
      mode = mode,
      barplot_base64 = result$barplot_base64,
      cnetplot_base64 = result$cnetplot_base64,
      unmapped = result$unmapped,
      csv_path = result$csv_path,
      zip_file = basename(result$zip_file),
      output_dir = result$output_dir
          

      
    ))
    
  }, error = function(e) {
    cat("❌ ERROR: ", e$message, "\n")
    res$status <- 500
    return(list(success = FALSE, error = e$message))
  })
}

########################################################
#KEGGMAP 
########################################################
#* @param mode The mode: "genelist" or "genelist_fc"
#* @param species The species: "mouse" or "human"
#* @param file_path Path to the input gene file
#* @param output_dir Path to the output folder
#* @post /keggmap
function(req, res) {
  cat("==== /keggmap CALLED ====\n")
  
  parsed <- tryCatch({ jsonlite::fromJSON(req$postBody) }, error = function(e) NULL)
  if (is.null(parsed) || is.null(parsed$file_path) || is.null(parsed$mode) ||
      is.null(parsed$species)) {
    res$status <- 400
    return(list(success = FALSE, error = "Missing required parameters"))
  }
  
  mode <- parsed$mode
  species <- parsed$species
  file_path <- parsed$file_path
  output_dir <- './Keggmap_output'


  # run KEGG analysis
  tryCatch({
    if (mode == "genelist") {
      cat("==== genelist mode ====\n")
      result <- run_kegg_from_gene_list(file_path, species)
    } else if (mode == "genelist_fc") {
      cat("==== genelist_fc mode ====\n")
      result <- run_kegg_from_gene_list_fc(file_path, species)
    } else {
      stop("Invalid mode. Use 'genelist' or 'genelist_fc'.")
    }

    return(list(
      success = TRUE,
      mode = mode,
      barplot_base64 = result$barplot_base64,
      top10_csv = result$top10_csv,
      full_csv = result$full_csv,
      unmapped = result$unmapped,
      maps = result$maps,
      foldchange_maps = result$foldchange_maps,
      zip_file = result$zip_file,     # ✅ path to zip file
      output_dir = result$output_dir
    ))
    
  }, error = function(e) {
    res$status <- 500
    return(list(success = FALSE, error = e$message))
  })
}

