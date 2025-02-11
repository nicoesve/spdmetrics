FROM rocker/r-ver:4.1.2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

# Install R packages
RUN R -e "install.packages(c('Matrix', 'expm', 'purrr', 'R6', 'languageserver'))"
