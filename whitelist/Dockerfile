FROM australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/whitelist_filter:kccg-chip-1.0

USER root

RUN Rscript -e 'if (!require(data.table)) { install.packages("data.table") }'
RUN Rscript -e 'if (!require(tidyr)) { install.packages("tidyr") }'
RUN Rscript -e 'if (!require(stringr)) { install.packages("stringr") }'