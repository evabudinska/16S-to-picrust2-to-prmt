# Data generation (PICRUSt2 native)

```bash
# Conda route (first time)
mamba env create -f scripts/picrust2/environment.yml
mamba activate picrust2-2.5.3

# Run
bash scripts/picrust2/run_picrust2.sh \
  --asv_table inst/extdata/toy_asv_table.tsv \
  --rep_seqs inst/extdata/toy_rep_seqs.fasta \
  --outdir results/picrust2

