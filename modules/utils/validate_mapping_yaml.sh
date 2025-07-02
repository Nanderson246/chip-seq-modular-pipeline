#!/bin/bash

# ------------------------------------------------------------
# Purpose: Validate and preview Sample_Type rules from YAML schema
# Usage:   bash modules/utils/validate_mapping_yaml.sh templates/mapping_schema.yaml
# Requires: yq (v4+, https://github.com/mikefarah/yq)
# Author: Nancy Anderson (concept with ChatGPT)
# ------------------------------------------------------------

set -euo pipefail

# === Input YAML file ===
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

# === Validate YAML format ===
echo "[INFO] Validating YAML syntax for: $YAML_FILE"
if ! yq eval '.' "$YAML_FILE" >/dev/null 2>&1; then
    echo "[ERROR] Invalid YAML syntax — check formatting!"
    exit 1
fi

# === Check for required sample_type_rules keys ===
echo "[INFO] Checking required keys under sample_type_rules..."

missing_keys=()
for key in ip_types control_types; do
    if ! yq eval ".sample_type_rules.$key" "$YAML_FILE" >/dev/null 2>&1; then
        missing_keys+=("$key")
    fi
done

if [[ ${#missing_keys[@]} -gt 0 ]]; then
    echo "[ERROR] Missing required keys under 'sample_type_rules': ${missing_keys[*]}"
    exit 1
fi

# === Extract and validate values ===
IP_TYPES=($(yq eval '.sample_type_rules.ip_types[]' "$YAML_FILE"))
CONTROL_TYPES=($(yq eval '.sample_type_rules.control_types[]' "$YAML_FILE"))

if [[ ${#IP_TYPES[@]} -eq 0 ]]; then
    echo "[ERROR] No entries found for .sample_type_rules.ip_types"
    exit 1
fi

if [[ ${#CONTROL_TYPES[@]} -eq 0 ]]; then
    echo "[ERROR] No entries found for .sample_type_rules.control_types"
    exit 1
fi

# === Report ===
echo "---------------------------------------------"
echo "✅ Allowed IP types:       ${IP_TYPES[*]}"
echo "✅ Allowed Control types:  ${CONTROL_TYPES[*]}"
echo "[SUCCESS] Sample type rules validated ✔️"
echo "---------------------------------------------"

