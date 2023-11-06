FROM debian:bookworm-slim

RUN apt-get update \
  && apt-get install -y perl \
  && rm -r /var/lib/apt/lists/* \
  && rm -r /var/cache/apt/* \
  && mkdir -p /annovar

COPY annovar/*.pl /annovar

ENV PATH=/annovar:${PATH}

WORKDIR /annovar