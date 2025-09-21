# R 4.5.1 (pairs with Bioconductor 3.21)
FROM rocker/r-ver:4.5.1
SHELL ["/bin/bash","-o","pipefail","-c"]
ENV CRAN=https://cloud.r-project.org

# ---- System deps ----
# - httpuv needs zlib, plumber needs libsodium
# - font libs help plotting (wordcloud, ggplot2, circlize)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential pkg-config \
    zlib1g-dev libsodium-dev \
    libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
    libpng-dev libjpeg-dev libtiff5-dev \
    libfreetype6-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    fonts-dejavu-core \
    curl \          
 && rm -rf /var/lib/apt/lists/*

# ---- Core CRAN helpers + plumber (install low-level deps explicitly) ----
RUN R -q -e "options(repos=c(CRAN=Sys.getenv('CRAN'))); \
             install.packages(c('remotes','httpuv','sodium','plumber'))"

# ---- Bioconductor 3.21 stack ----
RUN R -q -e "if (!requireNamespace('BiocManager', quietly=TRUE)) \
               install.packages('BiocManager', repos=Sys.getenv('CRAN')); \
             BiocManager::install(version='3.21', ask=FALSE, update=TRUE); \
             BiocManager::install(c( \
               'clusterProfiler','enrichplot','DOSE','fgsea', \
               'AnnotationDbi','GO.db','org.Hs.eg.db','org.Mm.eg.db', \
               'KEGGREST','GOSemSim','BiocParallel','pathview' \
             ), ask=FALSE, update=TRUE, dependencies=TRUE)"

# ---- CRAN packages used by your app ----
# (includes everything from your library() list: stringr, circlize, wordcloud, etc.)
RUN R -q -e "install.packages(c( \
  'jsonlite','dplyr','tidyr','readr','tibble','stringr','httr','xml2', \
  'openxlsx','zip','rentrez','base64enc','futile.logger', \
  'ggplot2','ggrepel','gridExtra','VennDiagram','UpSetR','wordcloud','RColorBrewer', \
  'purrr','circlize' \
), repos=Sys.getenv('CRAN'))"

# ---- Fail fast: verify that the important packages are installed & loadable ----
RUN R -q -e "req <- c( \
  # CRAN
  'plumber','jsonlite','dplyr','tidyr','readr','tibble','stringr','httr','xml2', \
  'openxlsx','zip','rentrez','base64enc','futile.logger', \
  'ggplot2','ggrepel','gridExtra','VennDiagram','UpSetR','wordcloud','RColorBrewer', \
  'purrr','circlize', \
  # Bioc
  'clusterProfiler','enrichplot','DOSE','fgsea','pathview', \
  'AnnotationDbi','GO.db','org.Hs.eg.db','org.Mm.eg.db','KEGGREST','GOSemSim','BiocParallel' \
); \
miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]; \
if (length(miss)) stop('Missing: ', paste(miss, collapse=', ')); \
invisible(lapply(req, function(p) library(p, character.only=TRUE))); \
cat('OK -> ', paste(req, sapply(req, function(p) as.character(packageVersion(p))), collapse=' | '), '\\n')"

# ---- App ----
WORKDIR /app
COPY . /app
EXPOSE 8000

# Run your Plumber API via your script
CMD ["Rscript","run_api.R"]
