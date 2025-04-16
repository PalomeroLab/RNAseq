#!/bin/bash
## bootstrap script for getting a fresh ubuntu instance up and running

# don't prompt for user input
export DEBIAN_FRONTEND=noninteractive
export NEEDRESTART_MODE=a

# update and install pre-requisites
sudo apt-get update && sudo apt-get upgrade -y
# sudo apt-get install -y unzip bzip2

# # install aws cli
# curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
# unzip awscliv2.zip && sudo ./aws/install && rm -rf ./awscliv2.zip && rm -rf ./aws
#
# # add completion
# echo "complete -C '/usr/local/bin/aws_completer' aws" >> ~/.bashrc

# Add Docker's official GPG key:
# sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
	"deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "${UBUNTU_CODENAME:-$VERSION_CODENAME}") stable" \
	| sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# configure no-sudo docker
sudo usermod -aG docker "$USER"
newgrp docker

# install micromamba
# curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
# micromamba shell init -s bash -r ~/.micromamba

# download the bioc run script
curl -fsSL https://raw.githubusercontent.com/Bioconductor/bioc-run/refs/heads/devel/bioc-run | sudo tee /usr/local/bin/bioc-run > /dev/null
sudo chmod +x /usr/local/bin/bioc-run

# restart shell
exec bash
