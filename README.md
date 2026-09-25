# merfish-mouse-testis

Nextflow pipeline for analysis of Vizgen MERSCOPE (MERFISH) data from mouse testis. Below is
a link to the report produced by this repo:

- **[Cell type annotation — b2r0_cellpose3d](reports/celltype_report_b2r0_report_gfm/b2r0_cellpose3d/celltype_report.md)**
  — QC, clustering, per-cell cell type calls, what they compose to per cluster, and the calls
  on tissue, for the lab's own cellpose segmentation of the B2 region.

## Analysis Overview

The pipeline reads a MERSCOPE region directory (or a cellpose re-segmentation of one) into a
SpatialData Zarr store (`create_spatialdata` / `create_spatialdata_cellpose`), clusters it on
the GPU across a Leiden resolution sweep (`cluster_spatialdata_gpu`), then correlates every
cell against a reference centroid table to call cell types (`annotate_celltypes`). Genes are
matched case-insensitively, so this mouse panel is annotated against a human reference atlas.
`create_centroids` builds per-cluster centroids for reporting, and a Quarto notebook renders
the cohort or per-sample results (`render_cohort` / `render_sample`).

Full per-step reference — inputs, outputs, and the reasoning behind each step's design — is
in `CLAUDE.md`.

## Running The Analysis

### 1. Clone the repository

```bash
git clone https://github.com/brent-biddy/merfish-mouse-testis.git
cd merfish-mouse-testis
```

### 2. Install Nextflow and Apptainer

Runs use a shared container via Apptainer, so no local Python environment is needed to run
the pipeline itself:

```bash
nextflow -version    # >=23.0
apptainer --version
```

See [these instructions](https://www.nextflow.io/docs/latest/install.html) for Nextflow and
[these instructions](https://apptainer.org/docs/admin/main/installation.html) for Apptainer
if either is missing.

### 3. Run the pipeline

Two entry points: `main.nf` chains the Vizgen path end to end in one invocation, `steps.nf`
runs one step at a time via `--step`.

```bash
# Full pipeline in one go
nextflow run main.nf -profile wsl --samplesheet assets/samplesheet.csv

# Or one step at a time — each artifact-producing step publishes a
# <step>_samplesheet.csv handoff sheet into outdir; point the next step at it.
nextflow run steps.nf --step create_spatialdata -profile wsl --samplesheet assets/samplesheet.csv
nextflow run steps.nf --step cluster_spatialdata_gpu -profile wsl \
    --samplesheet <results>/create_spatialdata_samplesheet.csv

# Verify wiring without executing anything
nextflow run steps.nf --step create_spatialdata -stub --samplesheet assets/samplesheet.csv
```

`-profile oscer` is the real target; `-profile wsl` is sized for a laptop or WSL2 box.
Outputs land in `results/<run_id>/<sample>/<step>/` — gitignored, beside the code, on both
profiles (only the work dir moves to scratch on `oscer`).

Full command reference, including the cellpose steps and the render steps, is in `CLAUDE.md`.

## Layout

```
├── steps.nf           single-step entry point; dispatches on --step
├── main.nf            full-pipeline entry point; chains create_spatialdata ->
│                      cluster_spatialdata_gpu -> annotate_celltypes -> create_centroids
├── nextflow.config    params, the defaults that always apply, and which site profiles exist
├── conf/              wsl.config / oscer.config — each states only its difference from the defaults
├── modules/           one .nf module per step, plus samplesheet.nf (the shared samplesheet reader)
├── bin/               python scripts invoked by the modules (argparse CLIs)
├── notebooks/         quarto notebooks rendered by the pipeline
├── assets/            sample sheets, the reference centroid tables, and the pptx template and
│                      lua filter a render needs
├── reports/           committed gfm renders, readable on GitHub (see the report above)
├── environment.yml    conda env for ad hoc, outside-the-pipeline analysis
└── data/raw/          raw instrument output (not committed)
```

`reports/` is the only output that is committed. Everything else under the repo is
gitignored and reproducible from `bin/` + `assets/` + `data/raw/`.
