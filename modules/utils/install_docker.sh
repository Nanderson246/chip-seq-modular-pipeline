#!/bin/bash
set -euo pipefail

echo "[INFO] Removing old Docker if present..."
sudo apt-get remove -y docker docker-engine docker.io containerd runc || true

echo "[INFO] Updating package index..."
sudo apt-get update -y

echo "[INFO] Installing required packages..."
sudo apt-get install -y \
    ca-certificates \
    curl \
    gnupg \
    lsb-release

echo "[INFO] Adding Docker's GPG key..."
sudo install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/$(. /etc/os-release && echo "$ID")/gpg \
  | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg

echo "[INFO] Setting up Docker's APT repository..."
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] \
  https://download.docker.com/linux/$(. /etc/os-release && echo "$ID") \
  $(lsb_release -cs) stable" \
  | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

echo "[INFO] Updating APT sources..."
sudo apt-get update -y

echo "[INFO] Installing latest Docker CE, CLI, and containerd..."
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

echo "[INFO] Enabling and starting Docker..."
sudo systemctl enable docker
sudo systemctl start docker

echo "[✅ DONE] Docker was successfully installed!"
docker --version
