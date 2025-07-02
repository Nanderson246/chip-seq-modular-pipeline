📁 Metadata Folder

This folder contains the metadata mapping file, which is critical for running the pipeline correctly.

📄 mapping.tsv
This is the main metadata file (currently empty).

You can find a template version in the templates/ folder.

Refer to the user guide for detailed instructions on how to fill in this file properly.

⚠️ Important Notes
This file is required for the pipeline to run. It defines how your samples are grouped, labeled, and processed.

It works in conjunction with:

The YAML configuration files in the templates/ folder

The FASTQ files in the samples/ folder

Together, these three components (mapping.tsv, *.yaml, and FASTQ files) are the foundation of your ChIP-seq or ATAC-seq analysis.
