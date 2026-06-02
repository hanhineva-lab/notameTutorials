FROM ghcr.io/bioconductor/bioconductor_docker:RELEASE_3_23-R-4.6.0

# Install system packages
RUN apt-get update && \
    apt-get -y install libcurl4-openssl-dev --no-install-recommends

# Initialize application project directory
WORKDIR /project
COPY renv.lock renv.lock

# Copy renv infrastructure
RUN mkdir -p renv
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# Restore R project library
RUN R -s -e "renv::restore()"

# Copy application files into image
COPY . .

