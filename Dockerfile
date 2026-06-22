# BioMate Docker Image
# Based on singularity.def configuration

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
    && rm -rf /var/lib/apt/lists/* \
    && echo "export BUILD_TIME=\"${BUILD_TIME}\"" > /etc/profile.d/build-time.sh

# Install build package
RUN python3 -m pip install --break-system-packages build

# Install BioMate package
COPY . /app
WORKDIR /app
RUN python3 -m pip install --break-system-packages .

# Set the default command
ENTRYPOINT ["biomate"]
CMD ["--help"]
