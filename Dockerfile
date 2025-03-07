FROM python:3.11-slim-bullseye AS base

ENV DEBIAN_FRONTEND=noninteractive

ARG BCFTOOLS_VERSION=${BCFTOOLS_VERSION:-1.21}
ARG ECHTVAR_VERSION=${ECHTVAR_VERSION:-v0.2.1}

RUN apt update && apt install -y --no-install-recommends \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        gnupg \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        openjdk-17-jdk-headless \
        wget \
        zip \
        zlib1g && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    pip install --no-cache-dir --upgrade pip

FROM base AS bcftools_compiler

RUN apt-get update && apt-get install --no-install-recommends -y \
        gcc \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        make \
        zlib1g-dev && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -xf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    strip bcftools plugins/*.so && \
    make DESTDIR=/bcftools_install install

FROM base AS base_bcftools

COPY --from=bcftools_compiler /bcftools_install/usr/local/bin/* /usr/local/bin/
COPY --from=bcftools_compiler /bcftools_install/usr/local/libexec/bcftools/* /usr/local/libexec/bcftools/

FROM base_bcftools AS base_bcftools_echtvar

ADD "https://github.com/brentp/echtvar/releases/download/${ECHTVAR_VERSION}/echtvar" /bin/echtvar

RUN chmod +x /bin/echtvar

COPY echtvar_config.json /echtvar_config.json

ENV ECHTVAR_CONFIG="/echtvar_config.json"

FROM base_bcftools_echtvar AS talos_gcloud

# Google Cloud SDK: use the script-based installation, as the Debian package is outdated.
ADD https://sdk.cloud.google.com install_glcoud.sh
RUN bash install_glcoud.sh --disable-prompts --install-dir=/opt && \
    rm install_glcoud.sh

ENV PATH=$PATH:/opt/google-cloud-sdk/bin

# Add in the additional requirements that are most likely to change.
COPY requirements*.txt README.md setup.py ./
RUN pip install --no-cache-dir -r requirements.txt
COPY src src/
RUN pip install --no-cache-dir --upgrade pip && pip install --no-cache-dir .[cpg]

FROM base_bcftools_echtvar AS talos_none

RUN echo "Skipping cloud dependency installation"

# Add in the additional requirements that are most likely to change.
COPY requirements*.txt README.md setup.py ./
RUN pip install --no-cache-dir -r requirements.txt
COPY src src/
RUN pip install --no-cache-dir .
