

library(plumber)
cat("🚀 Starting Plumber API\n")
dir.create("www", showWarnings = FALSE, recursive = TRUE)

# make sure we run from the script dir
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path)) setwd(dirname(normalizePath(script_path)))

pr <- plumber::plumb("plumber.R")
pr$run(host = "0.0.0.0", port = 8000)
