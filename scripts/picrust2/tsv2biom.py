
#!/usr/bin/env python3
"""
Convert a wide TSV (feature_id + sample columns) to BIOM (HDF5).
Assumes first column is feature/OTU/ASV ID, header row has sample IDs.
"""
import argparse
import pandas as pd
from biom.table import Table
from biom.util import biom_open  # <-- proper import

p = argparse.ArgumentParser()
p.add_argument("--input", required=True)
p.add_argument("--output", required=True)
args = p.parse_args()

df = pd.read_csv(args.input, sep="\t")
feat_col = df.columns[0]
df = df.set_index(feat_col)

# ensure numeric matrix
df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

table = Table(
    df.values,
    observation_ids=df.index.astype(str).tolist(),
    sample_ids=df.columns.astype(str).tolist()
)

with biom_open(args.output, 'w') as f:
    table.to_hdf5(f, "picrust2_to_prmt")

print(f"Wrote BIOM: {args.output}")


