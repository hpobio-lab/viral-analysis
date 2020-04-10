from rocker/tidyverse
WORKDIR /app

RUN Rscript -e 'install.packages("cowplot")'
RUN Rscript -e 'install.packages("GetoptLong")'
COPY scripts/variable_sites.R /usr/bin/
