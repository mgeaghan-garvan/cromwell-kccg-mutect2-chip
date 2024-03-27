FROM debian:bookworm-slim

ENV CONDA_ROOT_PREFIX /miniforge
ENV PATH $CONDA_ROOT_PREFIX/bin:/software:$PATH
ARG RVERSION=${RVERSION:-4.3.1}

RUN apt-get update && apt-get install -y libbz2-dev liblzma-dev git wget bash bzip2 zip curl gcc zlib1g-dev \
  && rm -r /var/lib/apt/lists/* \
  && rm -r /var/cache/apt/* \
  && mkdir -p /tmp/miniforge \
  && curl -L -o /tmp/miniforge/Miniforge3-Linux-x86_64.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
  && bash /tmp/miniforge/Miniforge3-Linux-x86_64.sh -b -p ${CONDA_ROOT_PREFIX} \
  && rm -rf /tmp/miniforge \
  && conda install -y --prefix ${CONDA_ROOT_PREFIX} -c conda-forge \
      r-base=${RVERSION} r-tidyverse r-optparse \
  && conda init

COPY annotate_chip annotate_chip_cohort /software/
RUN chmod +x /software/annotate_chip /software/annotate_chip_cohort
COPY *.R /R/scripts/