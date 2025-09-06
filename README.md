# 16S → PICRUSt2 → PRMT scoring (reproducible pipeline)

This repository demonstrates a fully reproducible workflow:

1. **16S amplicon data** → PICRUSt2 functional predictions 
2. **EC/KO → PRMT mapping** (preparation step) 
3. **PICRUSt2 outputs + mapping** → PRMT scores (R pipeline) 


## Dependencies

To run PICRUSt2 locally you need:

- **Conda (or Mamba)**: recommended package manager to create the PICRUSt2 environment  
  - Install Miniconda: https://docs.conda.io/en/latest/miniconda.html  
  - or Mambaforge: https://github.com/conda-forge/miniforge#mambaforge

- **R ≥ 4.2** with packages listed in `R/01_packages.R`

---

### Setup PICRUSt2 environment

Once conda/mamba is installed:

```bash
# create environment from spec
conda env create -f scripts/picrust2/environment.yml

# activate it
conda activate picrust2-2.5.3

# test biom + scikit-bio are available
python -c "from biom.table import Table; import skbio; print('biom + scikit-bio OK')"


## Step 1. 16S → PICRUSt2 
Use the wrapper script under [`scripts/picrust2/`](scripts/picrust2/) to generate PICRUSt2 outputs.

Example with the included toy dataset:


```bash
bash scripts/picrust2/run_picrust2.sh \
  --asv_table inst/extdata/toy_asv_table.tsv \
  --rep_seqs  inst/extdata/toy_rep_seqs.ids_are_sequences.fasta \
  --outdir    results/picrust2 \
  --threads 8 \
  --min_align 0.8 \
  --keep_tmp

```

  Outputs:

- results/picrust2/ec_pred_metagenome.tsv
- results/picrust2/ko_pred_metagenome.tsv
- results/picrust2/pathways_out/path_abun_unstrat.tsv
  
  
##   Step 2. EC/KO → PRMT mapping (preparation) 
PRMT mapping specifies how EC/KO features are aggregated into PRMTs.

A minimal example is provided: inst/extdata/prmt_mapping.tsv

To regenerate from KEGG (requires internet + KEGGREST):  


```bash
Rscript pipeline/build_prmt_mapping.R pipeline/config.yaml
```
  
This uses the kegg: section in pipeline/config.yaml and produces reaction_mapformula.tsv plus a PRMT mapping table.
 
 
##   Step 3. PICRUSt2 → PRMT 

The R pipeline consumes the PICRUSt2 outputs, plus metadata and mapping, and produces PRMT scores.

Quick start:

```bash
  Rscript pipeline/run_pipeline.R pipeline/config.yaml
```

Inputs (config):

- abundance_table: PICRUSt2 output (EC/KO/pathways table)
- metadata: sample metadata (must include sample_id)
- prmt_mapping: mapping of features → PRMTs

Outputs:

- results/prmt_scores.tsv, results/prmt_scores_z.tsv
- results/prmt_density.png, results/prmt_heatmap.png
- results/sessionInfo.txt



##  Configuration cheatsheet  

All pipeline parameters are centralized in pipeline/config.yaml.
Below is a commented example adapted to the toy dataset shipped in inst/extdata/.


# PICRUSt2 → PRMT pipeline config

inputs:

   PICRUSt2 unstratified table you want to use:
   
     - EC:  results/picrust2/ec_pred_metagenome.tsv
     - KO:  results/picrust2/ko_pred_metagenome.tsv
     - PWY: results/picrust2/pathways_out/path_abun_unstrat.tsv
  abundance_table: "results/picrust2/ec_pred_metagenome.tsv"

  # First column name in that table:
     EC/KO → "function"
     Pathways → "pathway"
  id_col: "function"

  # Sample metadata (must contain `sample_id`)
  metadata: "inst/extdata/metadata.tsv"

  # Mapping of PICRUSt2 features → PRMTs
  prmt_mapping: "inst/extdata/prmt_mapping.tsv"

processing:
  normalize: "relative"      # "relative" | "none"

scoring:
  method: "sum"              # "sum" | "mean"

output:
  dir: "results"             # all outputs go here
  prefix: "prmt"             # filename prefix
  mapformula_file: "reaction_mapformula.tsv"  # file created by mapping builder

kegg:
  kgml_dir: "kgml_files"      # where KGML XML files are cached
  scan_reactions: false       # true = crawl KEGG REST (slow, not needed normally)
  checkpoint: "processed_reactions.RDS"
  kos_to_rxns_file: "inst/extdata/reference/kos_rxns.tsv"  # KO↔Reaction map
  kolist: []                  # subset of KOs, e.g. ["K00001","K00002"]; empty = all
  all_kegg_prefix: "all_kegg" # prefix for *_KeggReactions.rda files



##  Toy run (EC → PRMT)  


This repo ships tiny toy inputs in inst/extdata/. To smoke-test the pipeline:


# (Optional) run PICRUSt2 native on the toy ASV/FASTA to produce outputs

```bash
bash scripts/picrust2/run_picrust2.sh \
  --asv_table inst/extdata/toy_asv_table.tsv \
  --rep_seqs  inst/extdata/toy_rep_seqs.fasta \
  --outdir    results/picrust2
```

# Configure pipeline to use EC predictions and toy metadata/mapping

```bash
cat > pipeline/config.yaml <<'YAML'
```

inputs:
  abundance_table: "results/picrust2/ec_pred_metagenome.tsv"
  id_col: "function"
  metadata: "inst/extdata/metadata.tsv"
  prmt_mapping: "inst/extdata/prmt_mapping.tsv"

processing:
  normalize: "relative"

scoring:
  method: "sum"

output:
  dir: "results"
  prefix: "prmt"
YAML

# Run R pipeline end-to-end

```bash
Rscript pipeline/run_pipeline.R pipeline/config.yaml
```



##################### Minimal workflow ##################### 

This repo is designed to be reproducible from raw ASV counts to PRMT scores.
Here is the shortest path to run everything end-to-end on the included toy data.


#### 1. Generate PICRUSt2 predictions

Run the wrapper script (inside Conda env or Docker/Singularity, see scripts/picrust2/README.md

```bash
bash scripts/picrust2/run_picrust2.sh \
  --asv_table inst/extdata/toy_asv_table.tsv \
  --rep_seqs  inst/extdata/toy_rep_seqs.fasta \
  --outdir    results/picrust2
```

This creates:

- results/picrust2/ec_pred_metagenome.tsv
- results/picrust2/ko_pred_metagenome.tsv
- results/picrust2/pathways_out/path_abun_unstrat.tsv

## 2. Adjust configuration

Edit pipeline/config.yaml so that:

##
inputs:
  abundance_table: "results/picrust2/ec_pred_metagenome.tsv"
  id_col: "function"
  metadata: "inst/extdata/metadata.tsv"
  prmt_mapping: "inst/extdata/prmt_mapping.tsv"

##
Other settings (normalization, scoring) can be left as defaults.

## 3. Run PRMT scoring pipeline

```bash
Rscript pipeline/run_pipeline.R pipeline/config.yaml
```

Outputs will appear in results/:

- prmt_scores.tsv
- prmt_scores_z.tsv
- prmt_density.png, prmt_heatmap.png
- sessionInfo.txt

#### 4. (Optional) Build KEGG-based mapping

If you want to regenerate the PRMT mapping from KEGG data:

```bash
Rscript pipeline/build_prmt_mapping.R pipeline/config.yaml
```

#### 5. Run tests

To verify the installation works:

```bash
Rscript tests/run_tests.R
```

This runs a small testthat suite on the toy dataset.


## Example data

We ship a minimal toy dataset under `inst/extdata/`:

| File                   | Description                                      |
|-------------------------|--------------------------------------------------|
| `toy_asv_table.tsv`     | Toy ASV abundance table (features × samples)     |
| `toy_rep_seqs.fasta`    | Representative 16S sequences for toy ASVs        |
| `metadata.tsv`          | Minimal sample metadata (2 samples, 1 grouping)  |
| `prmt_mapping.tsv`      | Example mapping of ECs to PRMTs                  |

Additionally, `scripts/picrust2/params.toy.json` contains an example JSON config for the PICRUSt2 wrapper script.




# After this, you can run:

bash scripts/picrust2/run_picrust2.sh \
  --asv_table inst/extdata/toy_asv_table.tsv \
  --rep_seqs  inst/extdata/toy_rep_seqs.fasta \
  --outdir    results/picrust2
```

Alternative: container

If you cannot/don’t want to install conda locally, use the Docker/Singularity mode:

```bash
bash scripts/picrust2/run_picrust2.sh \
  --asv_table inst/extdata/toy_asv_table.tsv \
  --rep_seqs  inst/extdata/toy_rep_seqs.fasta \
  --outdir    results/picrust2 \
  --container docker
```

