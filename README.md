# 🧬 ChIP-seq Modular Analysis Pipeline

This project provides a **modular, reproducible ChIP-seq pipeline** designed for human and mouse samples. It supports both **MACS3** and **HOMER** peak callers and offers a streamlined two-stage workflow:

1. **Pipeline 1** – Preprocessing and BAM cleaning  
2. **Pipeline 2** – Peak calling, reproducibility (IDR), and annotation

Inspired by the ENCODE standard, it is aimed at advanced users who want flexibility and control. The pipeline is written in **Bash**, with optional **Docker** support for consistent environments.

> ⚠️ Note: macOS is not supported. This project is tested only on Linux with POSIX-compliant tools.

---

## 📦 Features

- **Adapter trimming**, alignment with BWA, and BAM cleanup
- **Replicate QC** using deepTools
- **Peak calling** with MACS3 or HOMER
- **Signal generation** (fold change, p-value)
- **IDR analysis** with patched NumPy-compatible `idr`
- **Modular execution**: run full pipelines or individual steps
- **Metadata-driven**: based on a user-provided `mapping.tsv` file

---

## 🚀 Quick Start

1. Clone the repository  
2. Prepare your FASTQ samples in the `samples/` folder  
3. Fill out the `metadata/mapping.tsv` file  
4. Set up the genome and spike-in references  
5. Run:

```bash
bash run_pipeline1.sh
bash run_pipeline2.sh --caller macs3 --genome hg38

> 📘 **Looking for full documentation?**  
> Visit the [Wiki](https://github.com/Nanderson246/chip-seq-modular-pipeline/wiki) for setup, usage, and module details.

