#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 --asv_table <tsv> --rep_seqs <fasta> --outdir <dir> [--refdb <dir>] [--container docker|singularity] [--threads N] [--min_align 0.8] [--keep_tmp]"
  exit 1
}

ASV_TABLE=""
REP_SEQS=""
OUTDIR="results/picrust2"
REFDB=""
CONTAINER="none"
THREADS=2
MIN_ALIGN="0.8"
KEEP_TMP=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --asv_table) ASV_TABLE="$2"; shift 2;;
    --rep_seqs)  REP_SEQS="$2";  shift 2;;
    --outdir)    OUTDIR="$2";    shift 2;;
    --refdb)     REFDB="$2";     shift 2;;
    --container) CONTAINER="$2"; shift 2;;
    --threads)   THREADS="$2";   shift 2;;
    --min_align) MIN_ALIGN="$2"; shift 2;;
    --keep_tmp)  KEEP_TMP=1;     shift 1;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

[[ -f "$ASV_TABLE" ]] || { echo "ASV table not found: $ASV_TABLE"; exit 2; }
[[ -f "$REP_SEQS"  ]] || { echo "Rep seqs not found: $REP_SEQS"; exit 2; }
mkdir -p "$OUTDIR"

WORKDIR="$(mktemp -d)"
if [[ $KEEP_TMP -eq 0 ]]; then
  trap 'rm -rf "$WORKDIR"' EXIT
else
  echo "[INFO] Keeping temp dir: $WORKDIR"
fi

echo "[INFO] Temp dir: $WORKDIR"
echo "[INFO] Converting TSV -> BIOM ..."
python "$(dirname "$0")/tsv2biom.py" --input "$ASV_TABLE" --output "$WORKDIR/table.biom"
echo "Wrote BIOM: $WORKDIR/table.biom"

# build picrust2 command
PIC_CMD=(picrust2_pipeline.py -s "$REP_SEQS" -i "$WORKDIR/table.biom" -o "$WORKDIR/out" -p "$THREADS" --min_align "$MIN_ALIGN")
[[ -n "$REFDB" ]] && PIC_CMD+=(--custom_ref_dir "$REFDB")

run_local() { "${PIC_CMD[@]}"; }
run_docker() {
  local IMG="biocontainers/picrust2:v2.5.3_cv1"
  docker run --rm -u $(id -u):$(id -g) \
    -v "$PWD":"$PWD" -v "$WORKDIR":"$WORKDIR" -w "$PWD" \
    "$IMG" "${PIC_CMD[@]}"
}
run_sing() {
  local IMG="biocontainers/picrust2:v2.5.3_cv1"
  singularity exec docker://"$IMG" "${PIC_CMD[@]}"
}

echo "[INFO] Running PICRUSt2 (threads=$THREADS, min_align=$MIN_ALIGN, container=$CONTAINER) ..."
case "$CONTAINER" in
  none)        run_local;;
  docker)      run_docker;;
  singularity) run_sing;;
  *)           echo "Unknown --container $CONTAINER"; exit 3;;
esac

# quick listing for debug
echo "[INFO] PICRUSt2 output tree (depth 2):"
find "$WORKDIR/out" -maxdepth 2 -type f | sed "s|$PWD/||"

# helper to copy .tsv or .tsv.gz (gunzip to .tsv)
copy_tsv_like() {
  local src_dir="$1" ; local base="$2" ; local dest="$3"
  if [[ -f "$src_dir/$base.tsv" ]]; then
    cp "$src_dir/$base.tsv" "$dest"
    return 0
  elif [[ -f "$src_dir/$base.tsv.gz" ]]; then
    gunzip -c "$src_dir/$base.tsv.gz" > "$dest"
    return 0
  fi
  return 1
}

# collect outputs (decompress if needed)
ok=0
if copy_tsv_like "$WORKDIR/out/KO_metagenome_out" "pred_metagenome_unstrat" "$OUTDIR/ko_pred_metagenome.tsv"; then
  echo "[INFO] Wrote $OUTDIR/ko_pred_metagenome.tsv"; ok=1
fi
if copy_tsv_like "$WORKDIR/out/EC_metagenome_out" "pred_metagenome_unstrat" "$OUTDIR/ec_pred_metagenome.tsv"; then
  echo "[INFO] Wrote $OUTDIR/ec_pred_metagenome.tsv"; ok=1
fi
mkdir -p "$OUTDIR/pathways_out"
if copy_tsv_like "$WORKDIR/out/pathways_out" "path_abun_unstrat" "$OUTDIR/pathways_out/path_abun_unstrat.tsv"; then
  echo "[INFO] Wrote $OUTDIR/pathways_out/path_abun_unstrat.tsv"; ok=1
fi

# provenance
{
  echo "picrust2: $(picrust2_pipeline.py --version 2>/dev/null || echo 'unknown')"
  echo "date: $(date -u +%F'T'%H:%M:%SZ)"
  echo "threads: $THREADS"
  echo "min_align: $MIN_ALIGN"
} > "$OUTDIR/software_versions.yml"

jq -n \
  --arg asv "$ASV_TABLE" \
  --arg rep "$REP_SEQS" \
  --arg out "$OUTDIR" \
  --arg refdb "${REFDB:-}" \
  --arg threads "$THREADS" \
  --arg min_align "$MIN_ALIGN" \
  '{asv_table:$asv, rep_seqs:$rep, outdir:$out, refdb:$refdb, threads:$threads, min_align:$min_align}' \
  > "$OUTDIR/params_used.json"

if [[ $ok -eq 0 ]]; then
  echo "[WARN] Expected PICRUSt2 outputs not found. Check logs under: $WORKDIR/out/logs"
  echo "       You can re-run with --keep_tmp to inspect."
fi

echo "PICRUSt2 done → $OUTDIR"
BASH

chmod +x scripts/picrust2/run_picrust2.sh

