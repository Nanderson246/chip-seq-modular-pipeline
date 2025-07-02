#!/usr/bin/env python3

#use: python3 modules/utils/validate_mapping.py --mapping metadata/mapping.tsv --schema templates/mapping_schema.yaml


import argparse
import csv
import yaml
import sys
import os
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Validate a mapping file against a YAML schema.")
    parser.add_argument("--mapping", required=True, help="Path to the mapping file (.tsv, .csv, or .txt)")
    parser.add_argument("--schema", required=True, help="Path to the YAML schema file")
    return parser.parse_args()

def detect_delimiter(filename):
    ext = os.path.splitext(filename)[1].lower()
    if ext == ".csv":
        return ","
    return "\t"  # default to TSV or TXT

def load_schema(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)

def validate_header(header, schema):
    missing = [field for field in schema.get("required_fields", []) if field not in header]
    if missing:
        print(f"[ERROR] Missing required columns: {', '.join(missing)}")
        return False
    return True

def validate_rows(reader, header, schema):
    field_rules = schema.get("field_rules", {})
    valid = True
    for i, row in enumerate(reader, start=2):
        row_dict = dict(zip(header, [cell.strip() for cell in row]))
        for field, rules in field_rules.items():
            val = row_dict.get(field, "").strip()

            # Check allowed values
            if "allowed" in rules:
                if val and val not in rules["allowed"]:
                    print(f"[ERROR] Line {i}: Invalid value '{val}' in '{field}'. Allowed: {rules['allowed']}")
                    valid = False

            # Check regex pattern
            if "regex" in rules:
                pattern = rules["regex"]
                if val and not re.match(pattern, val):
                    print(f"[ERROR] Line {i}: Value '{val}' in '{field}' does not match pattern '{pattern}'")
                    valid = False
    return valid

def check_empty_lines(filename):
    with open(filename, "r") as f:
        for i, line in enumerate(f, start=1):
            if not line.strip():
                print(f"[WARNING] Empty line detected at line {i}")

def main():
    args = parse_args()
    check_empty_lines(args.mapping)  # Add this early
    delim = detect_delimiter(args.mapping)
    schema = load_schema(args.schema)

    with open(args.mapping, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=delim)
        try:
            header = next(reader)
            header = [h.strip() for h in header]
        except StopIteration:
            print("[ERROR] Mapping file is empty.")
            sys.exit(1)

        if not validate_header(header, schema):
            sys.exit(1)

        if not validate_rows(reader, header, schema):
            sys.exit(1)

    print("[INFO] ✅ Mapping file passed validation.")

if __name__ == "__main__":
    main()

