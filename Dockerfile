FROM python:3.11-bullseye AS base

COPY --from=ghcr.io/astral-sh/uv:0.5.5 /uv /uvx /bin/

# Set up working directory for the project
WORKDIR /talos

RUN apt update && apt install -y --no-install-recommends \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        git \
        gnupg \
        openjdk-11-jdk-headless \
        wget \
        zip && \
    apt clean

# Install the project's dependencies using the lockfile and settings
# Copy `pyproject.toml` and `uv.lock` into `/talos` explicitly
COPY pyproject.toml uv.lock /talos/

# Install the project's dependencies using the lockfile and settings
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --frozen --no-install-project

# Add the project source code from src/cpg-flow
ADD . /talos
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen

# Place executables in the environment at the front of the path
ENV PATH="/talos/.venv/bin:$PATH"
ENV PYTHONPATH="/talos:${PYTHONPATH}"

# install nextflow
ADD https://get.nextflow.io nextflow
RUN chmod +x nextflow && \
    mv nextflow /usr/bin && \
    nextflow self-update


FROM base AS talos_gcloud

# Google Cloud SDK: use the script-based installation, as the Debian package is outdated.
ADD https://sdk.cloud.google.com install_gcloud.sh
RUN bash install_gcloud.sh --disable-prompts --install-dir=/opt && \
    rm install_gcloud.sh

ENV PATH=$PATH:/opt/google-cloud-sdk/bin


FROM base AS talos_none

RUN echo "Skipping cloud dependency installation"
