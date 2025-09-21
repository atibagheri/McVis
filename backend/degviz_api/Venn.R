library(VennDiagram)
library(UpSetR)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(base64enc)


generate_venn_or_upset <- function(paths, output_png = "venn_output.png", output_pdf = "venn_output.pdf") {
  # ✅ Use the provided file paths directly
  files <- paths
  cat("✅ Files detected:\n"); print(files)
  
  if (length(files) < 2) {
    stop("❌ Please provide at least 2 TXT files.")
  }
  
  # ✅ Read datasets
  datasets <- lapply(files, function(f) {
    df <- read.table(f, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    if (!("Gene.ID" %in% colnames(df))) {
      stop(paste("❌ File", f, "does not contain a Gene.ID column"))
    }
    unique(df$Gene.ID)
  })
  
  names(datasets) <- basename(tools::file_path_sans_ext(files))
  num_sets <- length(datasets)
  cat("✅ Number of sets:", num_sets, "\n")
  
  # ✅ Plot logic
  plot_logic <- function() {
    if (num_sets == 2) {
      a1 <- length(datasets[[1]])
      a2 <- length(datasets[[2]])
      inter12 <- length(intersect(datasets[[1]], datasets[[2]]))
      draw.pairwise.venn(
        area1 = a1, area2 = a2, cross.area = inter12,
        category = names(datasets),
        fill = c("#56B4E9", "#E69F00"),
        cat.col = c("#56B4E9", "#E69F00"),
        cat.cex = 3, cex = 3, lty = "solid", scaled = FALSE,  
        cat.pos  = c(0, 0),        # 270° = left, 90° = right
        cat.dist = c(0.015, 0.015),     # push outward a bit
        cat.just = list(c(0.5,0),    # left label right-justified
                  c(0.5, 0))
      )
    } else if (num_sets == 3) {
      a1 <- length(datasets[[1]])
      a2 <- length(datasets[[2]])
      a3 <- length(datasets[[3]])
      n12 <- length(intersect(datasets[[1]], datasets[[2]]))
      n23 <- length(intersect(datasets[[2]], datasets[[3]]))
      n13 <- length(intersect(datasets[[1]], datasets[[3]]))
      n123 <- length(Reduce(intersect, datasets[1:3]))
      cat("Set sizes:", a1, a2, a3, "\n")
      cat("Pairwise intersections:", n12, n23, n13, "\n")
      cat("Triple intersection:", n123, "\n")
      
      draw.triple.venn(
        area1 = 100, area2 = 100, area3 = 100,
        n12 = n12, n23 = n23, n13 = n13, n123 = n123,
        category = names(datasets),
        fill = c("#56B4E9", "#E69F00", "#F0E442"),
        cat.col = c("#56B4E9", "#E69F00", "#F0E442"),
        cat.cex = 3, cex = 3, lty = "solid", 
        scaled = FALSE,
        cat.pos = c(-20, 20, 180),   # tweak until you like
        cat.dist = c(0.1, 0.1, 0.1)
      )
    } else if (num_sets == 4 || num_sets == 5) {
      fill_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")[1:num_sets]
      # Generate in-memory Venn object, no log file
      venn.plot <- venn.diagram(
        x = datasets,
        category.names = names(datasets),
        filename = NULL,     # don't create a file
        output = FALSE,      # suppress internal writing
        col = "blue",
        disable.logging = TRUE,
        fill = fill_colors,
        alpha = 0.5,
        cex = 3,
        cat.cex = 2,
        cat.dist = 0.09,
        scaled = FALSE,
        main = NULL,
        print.mode = "raw"
      )
      # Render directly to the PNG/PDF device
      grid.newpage()
      grid.draw(venn.plot)
     
    } else if (num_sets >= 6) {
      df_list <- lapply(names(datasets), function(name) {
        data.frame(Gene.ID = datasets[[name]], Set = name)
      })
      combined_df <- bind_rows(df_list)
      df_wide <- combined_df %>%
        mutate(value = 1) %>%
        tidyr::spread(Set, value, fill = 0)
      grid.newpage()
      print(upset(df_wide, sets = names(datasets), main.bar.color = "skyblue", point.size = 3, order.by = "freq"))
    } else {
      stop("❌ Something went wrong.")
    }
  }
  
  # ✅ Save to PNG
  png(output_png, width = 800, height = 800)
  plot_logic()
  dev.off()
  
  # ✅ Save to PDF
  pdf(output_pdf, width = 12, height = 8)
  plot_logic()
  dev.off()
  
  message("✅ Overlap plot saved as PNG: ", output_png)
  message("✅ Overlap plot saved as PDF: ", output_pdf)
  
  # ✅ Return base64
  img_bytes <- readBin(output_png, what = "raw", n = file.info(output_png)$size)
  img_base64 <- base64encode(img_bytes)
  img_base64 <- base64enc::base64encode(img_bytes)

  
  return(list(
    output_png = normalizePath(output_png),
    output_pdf = normalizePath(output_pdf),
    png_base64 = img_base64
  ))
}
