library(plumber)

cat("🚀 Starting Plumber API\n")
pr <- plumber::plumb("plumber.R")
# run on port 8000
pr$run(host = "0.0.0.0", port = 8000)
