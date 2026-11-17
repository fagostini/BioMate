# BioMate Docker Image
# Based on singularity.def configuration

# Build the optional Rust kernels (biomate_kernels: nspector and index --engine rust)
FROM ubuntu:24.04 AS kernels
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y curl ca-certificates build-essential python3 python3-venv \
    && rm -rf /var/lib/apt/lists/*
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal
ENV PATH="/root/.cargo/bin:$PATH"
RUN python3 -m venv /opt/maturin && /opt/maturin/bin/pip install "maturin>=1.5,<2.0"
COPY rust/biomate_kernels /src/biomate_kernels
WORKDIR /src/biomate_kernels
RUN /opt/maturin/bin/maturin build --release --out /wheels

FROM ubuntu:24.04

# Avoid interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Set build time as environment variable
RUN BUILD_TIME=$(date) && \
    apt-get update && \
    apt-get install -y \
        python3 \
        python3-pip \
        build-essential \
        python3-dev \
        python3-venv \
    && rm -rf /var/lib/apt/lists/* \
    && echo "export BUILD_TIME=\"${BUILD_TIME}\"" > /etc/profile.d/build-time.sh

# Create and activate Python virtual environment
RUN python3 -m venv /app/venv
ENV PATH="/app/venv/bin:$PATH"

# Install build tools in venv
RUN pip install --upgrade pip setuptools wheel
RUN pip install build

# Install BioMate package from local source
COPY . /app
WORKDIR /app
RUN pip install .
COPY --from=kernels /wheels /tmp/wheels
RUN pip install /tmp/wheels/*.whl && rm -r /tmp/wheels

# Web interface port (see docker-compose.yml)
EXPOSE 8765
