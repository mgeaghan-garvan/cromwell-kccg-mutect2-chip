FROM australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/whitelist_filter:kccg-chip-1.0

USER root

RUN Rscript -e 'if (!require(tidyr)) { install.packages("tidyr") }'
RUN Rscript -e 'if (!require(stringr)) { install.packages("stringr") }'
RUN cd /usr/local/bin && wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && mv bedtools.static.binary bedtools && chmod +x bedtools
