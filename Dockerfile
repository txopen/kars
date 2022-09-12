# Install base image
ARG VERSION="4.1.0"
FROM rocker/r-ver:$VERSION as base

# Install needed packages and git
RUN apt-get update
RUN apt-get -y install libgit2-dev 
RUN apt-get -y install libcurl4-openssl-dev 
RUN apt-get -y install libssl-dev
RUN apt-get -y install r-base

# Install R packages
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "devtools::install_github('https://github.com/txopen/histoc')"
RUN Rscript -e "devtools::install_version('shinyjs', version = '2.1.0', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "devtools::install_version('shinythemes', version = '1.2.0', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "devtools::install_version('shinycssloaders', version = '1.0.0', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "devtools::install_version('shinybusy', version = '0.3.1', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "devtools::install_version('synchronicity', version = '1.3.5', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "devtools::install_version('DT', version = '0.24', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "devtools::install_version('gtsummary', version = '1.6.1', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "devtools::install_version('gt', version = '0.7.0', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "devtools::install_version('tidyverse', version = '1.3.2', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "devtools::install_version('data.table', version = '1.14.2', repos = 'http://cran.us.r-project.org')"

# Copy files
RUN mkdir www
COPY www www
COPY server.R .
COPY ui.R .

EXPOSE 3838
# Run app
CMD ["Rscript", "-e", "shiny::runApp()"]