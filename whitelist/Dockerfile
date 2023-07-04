FROM australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/whitelist_filter:kccg-chip-1.0

USER root

RUN Rscript -e 'if (!require(tidyr)) { install.packages("tidyr") }'
RUN Rscript -e 'if (!require(stringr)) { install.packages("stringr") }'
RUN apt-get update && apt-get install -y libbz2-dev liblzma-dev \
  && cd /usr/local/bin \
  && wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary \
  && mv bedtools.static.binary bedtools \
  && chmod +x bedtools \
  && mkdir -p /software \
  && cd /software \
  && wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 \
  && wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 \
  && tar -xjvf htslib-1.17.tar.bz2 \
  && tar -xjvf bcftools-1.17.tar.bz2 \
  && cd /software/htslib-1.17 \
  && ./configure --prefix=/software/htslib-1.17/build \
  && make \
  && make install \
  && cd - \
  && cd /software/bcftools-1.17 \
  && ./configure --prefix=/software/bcftools-1.17/build \
  && make \
  && make install \
  && cd -
ENV PATH="${PATH}:/software/htslib-1.17/build/bin:/software/bcftools-1.17/build/bin"
