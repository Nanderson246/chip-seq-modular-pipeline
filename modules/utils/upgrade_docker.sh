#!/bin/bash
set -euo pipefail

echo "🔍 Checking for existing Docker installation..."
if command -v docker &>/dev/null; then
  echo "✅ Docker version: $(docker --version)"
else
  echo "⚠️  Docker not currently installed. This script will install the latest version."
fi

echo "🔄 Removing old Docker versions (if any)..."
sudo apt remove -y docker docker-engine docker.io containerd runc || true

echo "📦 Installing dependencies..."
sudo apt update
sudo apt install -y ca-certificates curl gnupg lsb-release

echo "🔑 Adding Docker’s official GPG key..."
sudo install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | \
  sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg

echo "📚 Setting up Docker stable repository..."
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] \
  https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

echo "🔄 Updating package list..."
sudo apt update

echo "⬇️ Installing latest Docker engine and CLI..."
sudo apt install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

echo "✅ Docker successfully updated to:"
docker --version

