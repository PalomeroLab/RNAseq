#!/bin/bash
## This script sets up a Docker container for RStudio Server with Bioconductor.
## Original script can be found at `<https://github.com/Bioconductor/bioc-run>`
#!/usr/bin/env bash

# Path to the RStudio home directory inside the container, mapped from host
DOCKER_HOME="$HOME/rstudio"

# Path to shared R package library between host and container
DOCKER_RPKGS="$HOME/.my-bioc-packages"

# Port on the host that maps to container's RStudio Server port (8787)
PORT=8787

# Password for RStudio login (user is 'rstudio')
PASSWORD=bioc

# Bioconductor Docker image to use
IMAGE="bioconductor/bioconductor_docker:RELEASE_3_20"

# Name to assign to the Docker container instance
CONTAINER_NAME="bioconductor-rstudio"

# Create mount directories if they don't exist
mkdir -p "$DOCKER_HOME" "$DOCKER_RPKGS"

# You only need to run this once to set up the container
# Since we have `--restart always`, it should start automatically on boot
docker run -d \
  --name "$CONTAINER_NAME" \
  -e PASSWORD="$PASSWORD" \
  -e USERID="$(id -u)" \
  -v "$DOCKER_HOME":/home/rstudio \
  -v "$DOCKER_RPKGS":/usr/local/lib/R/host-site-library \
  -v "$HOME/.local/bin":/home/rstudio/.local/bin \
  -p "$PORT":8787 \
  --restart always \
  "$IMAGE" \
  bash -c 'export PATH="/home/rstudio/.local/bin:$PATH" && /init'
