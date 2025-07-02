#!/bin/bash

# Purpose: Validate and read YAML metadata rules (for Sample_Type checking)
# Author: Nancy Anderson (concept by ChatGPT)
# Usage: bash validate_yaml.sh <your_yaml_file.yml>

set -euo pipefail

# === Input YAML ===
YAML_FILE="${1:-}"

if [[ -z "$YAML_FILE" ]]; then
    echo "[ERROR] No YAML file specified."
    echo "Usage: $0 <yaml_file>"
    exit 1
fi

if [[ ! -f "$YAML_FILE" ]]; then
    echo "[ERROR] YAML file not found: $YAML_FILE"
    exit 1
fi

# === Validate YAML structure ===
echo "[INFO] Validating YAML format..."
if ! yq eval '.' "$YAML_FILE" >/dev/null 2>&1; then
    echo "[ERROR] Invalid YAML format."
    exit 1
fi

# === Read allowed sample types ===
echo "[INFO] Extracting allowed Sample_Types..."

IP_TYPES=($(yq eval '.sample_type_rules.ip_types[]' "$YAML_FILE"))
CONTROL_TYPES=($(yq eval '.sample_type_rules.control_types[]' "$YAML_FILE"))

# === Print out what was read ===
echo "✅ Allowed IP types: ${IP_TYPES[*]}"
echo "✅ Allowed Control types: ${CONTROL_TYPES[*]}"

# === Optional: Check if lists are empty ===
if [[ ${#IP_TYPES[@]} -eq 0 || ${#CONTROL_TYPES[@]} -eq 0 ]]; then
    echo "[WARN] One or both allowed type lists are empty!"
else
    echo "[SUCCESS] YAML sample_type_rules properly loaded."
fi

