# Reproducing PICRUSt2 outputs (native)

This folder runs the official `picrust2_pipeline.py` on an **ASV table (TSV)** and **rep seqs (FASTA)** and collects the standard outputs used by the PRMT pipeline.

## Quick start (Conda)
```bash
# 1) Create env (once)
mamba env create -f environment.yml
mamba activate picrust2-2.5.3

# 2) Run on toy or real data
bash run_picrust2.sh \
  --asv_table inst/extdata/toy_asv_table.tsv \
  --rep_seqs inst/extdata/toy_rep_seqs.fasta \
  --outdir results/picrust2

