# Running PICRUSt2 via Docker / Singularity

This repository provides a native runner for PICRUSt2 (`scripts/picrust2/run_picrust2.sh`)
that can execute entirely inside a container — no local installs needed.

- **Container image:** `biocontainers/picrust2:v2.5.3_cv1` (version pinned)
- **Inputs:** TSV ASV table (wide; first column = feature_id), FASTA representative sequences
- **Outputs:** KO/EC/pathway tables compatible with the R PRMT pipeline

---

## Quick start (Docker)

```bash
# From the repo root
bash scripts/picrust2/run_picrust2.sh \
  --asv_table inst/extdata/toy_asv_table.tsv \
  --rep_seqs  inst/extdata/toy_rep_seqs.fasta \
  --outdir    results/picrust2 \
  --container docker

