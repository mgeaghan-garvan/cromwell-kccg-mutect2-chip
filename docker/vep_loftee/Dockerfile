FROM ensemblorg/ensembl-vep:release_103.1

USER root

RUN apt-get update && apt-get install -y wget libncurses5-dev
RUN mkdir -p /apps /plugins
RUN cd /apps && \
    wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
    tar -xjvf samtools-1.13.tar.bz2 && \
    cd samtools-1.13 && \
    mkdir -p bin && \
    ./configure --prefix=/apps/samtools-1.13/bin/ && \
    make && \
    make install && \
    cd ../ && \
    wget https://github.com/samtools/bcftools/releases/download/1.13/bcftools-1.13.tar.bz2 && \
    tar -xjvf bcftools-1.13.tar.bz2 && \
    cd bcftools-1.13 && \
    mkdir -p bin && \
    ./configure --prefix=/apps/bcftools-1.13/bin/ && \
    make && \
    make install && \
    cd ../ && \
    wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 && \
    tar -xjvf htslib-1.13.tar.bz2 && \
    cd htslib-1.13 && \
    mkdir -p bin && \
    ./configure --prefix=/apps/htslib-1.13/bin/ && \
    make && \
    make install && \
    cd ../
RUN cd /plugins && \
    wget https://github.com/konradjk/loftee/archive/refs/tags/v1.0.3.tar.gz && \
    tar -xzvf v1.0.3.tar.gz && \
    cp loftee-1.0.3/LoF.pm . && \
    cpanm DBD::SQLite

ENV PATH="/apps/samtools-1.13/bin/bin:/apps/bcftools-1.13/bin/bin:/apps/htslib-1.13/bin/bin:${PATH}"

USER vep