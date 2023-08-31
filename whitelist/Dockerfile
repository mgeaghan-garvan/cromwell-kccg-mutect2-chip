FROM debian:bookworm-slim

ENV MAMBA_ROOT_PREFIX /mambaforge
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH
ARG RVERSION=${RVERSION:-4.3.1}

RUN apt-get update && apt-get install -y libbz2-dev liblzma-dev git wget bash bzip2 zip curl gcc zlib1g-dev \
  && rm -r /var/lib/apt/lists/* \
  && rm -r /var/cache/apt/* \
  && mkdir -p /tmp/mambaforge \
  && curl -L -o /tmp/mambaforge/Mambaforge-Linux-x86_64.sh https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh \
  && bash /tmp/mambaforge/Mambaforge-Linux-x86_64.sh -b -p ${MAMBA_ROOT_PREFIX} \
  && rm -rf /tmp/mambaforge \
  && mamba install -y --prefix ${MAMBA_ROOT_PREFIX} -c conda-forge -c bioconda \
      r-base=${RVERSION} r-tidyverse \
  && mamba init \
  && cd /usr/local/bin \
  && wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary \
  && mv bedtools.static.binary bedtools \
  && chmod +x bedtools \
  && mkdir -p /software/tmp \
  && mkdir -p /software/htslib-1.17 \
  && mkdir -p /software/bcftools-1.17 \
  && cd /software/tmp \
  && wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 \
  && wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 \
  && tar -xjvf htslib-1.17.tar.bz2 \
  && tar -xjvf bcftools-1.17.tar.bz2 \
  && cd /software/tmp/htslib-1.17 \
  && ./configure --prefix=/software/htslib-1.17/build \
  && make \
  && make install \
  && cd - \
  && cd /software/tmp/bcftools-1.17 \
  && ./configure --prefix=/software/bcftools-1.17/build \
  && make \
  && make install \
  && cd /software \
  && rm -rf /software/tmp \
  && mkdir /whitelist

ENV PATH="${PATH}:/software/htslib-1.17/build/bin:/software/bcftools-1.17/build/bin"

COPY *.R import/ chip_variant_definitions.csv somaticism_filter_transcripts.txt /whitelist/

WORKDIR /whitelist