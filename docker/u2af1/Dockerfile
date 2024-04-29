FROM rust:slim-bookworm

RUN apt-get update \
  && apt-get install -y git build-essential gawk \
  && rm -r /var/lib/apt/lists/* \
  && rm -r /var/cache/apt/* \
  && mkdir /software \
  && cd /software \
  && git clone https://github.com/weinstockj/pileup_region.git \
  && cd pileup_region \
  && cargo build --release \
  && ln -s /software/pileup_region/target/release/pileup_region /usr/local/bin/pileup_region