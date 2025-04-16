#!/bin/bash
## This script installs micromamba and the AWS CLI on a fresh Ubuntu instance.
## No superuser privileges are required.

set -euo pipefail

BIN_DIR="$HOME/.local/bin"
MAMBA_ROOT="$HOME/.micromamba"
BASHRC="$HOME/.bashrc"

mkdir -p "$BIN_DIR"

# Install micromamba
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
	| tar -xvj -C "$BIN_DIR" --strip-components=1 bin/micromamba

# Init micromamba shell
"$BIN_DIR/micromamba" shell init -s bash -r "$MAMBA_ROOT"

# Install AWS CLI
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip -q awscliv2.zip
./aws/install -i "$HOME/.local/aws-cli" -b "$BIN_DIR"
rm -rf aws awscliv2.zip

# Ensure PATH and autocompletion are in .bashrc
grep -qxF "export PATH=\"$BIN_DIR:\$PATH\"" "$BASHRC" \
	|| echo "export PATH=\"$BIN_DIR:\$PATH\"" >> "$BASHRC"

grep -qxF "complete -C '$BIN_DIR/aws_completer' aws" "$BASHRC" \
	|| echo "complete -C '$BIN_DIR/aws_completer' aws" >> "$BASHRC"

echo "Done. Restarting shell..."

exec bash
