# merfish-mouse-testis

Analysis of Vizgen MERSCOPE (MERFISH) data from mouse testis, annotated against a human
reference atlas (genes matched case-insensitively across species).

## Setup

There is no project-specific environment. Everything in `bin/` runs in a shared container,
and each `conf/<profile>.config` names two: a CPU image for every step, and a GPU image for the
one process that declares `label 'gpu'`.

```bash
apptainer exec oras://ghcr.io/brent-biddy/python_cpu-sif:1.0.2 bin/create_spatialdata.py --help

apptainer exec oras://ghcr.io/brent-biddy/python_gpu-sif:1.0.2 bin/cluster_spatialdata_gpu.py --help
```

Both are published as `oras://` artifacts, which are already SIFs, so a pull skips the
OCI-to-SIF conversion a `docker://` reference would pay for.

## Running the pipeline

Two entry scripts, and they interoperate.

```bash
# One step at a time, while the analysis is being worked out.
nextflow run steps.nf -profile wsl --step create_spatialdata --samplesheet assets/samplesheet.csv

# The Vizgen path chained end to end, as one run.
nextflow run main.nf -profile wsl --samplesheet assets/samplesheet.csv
```

Steps are selected by name; `steps.nf` chains nothing, so any step can be rerun on its own.
`main.nf` calls the same module workflows in order — an example of the shape rather than a
settled analysis. Because every step writes the next one's samplesheet into
`results/<run_id>/`, a chained run leaves the same breadcrumbs a stepwise one does: run the
whole thing, then re-run a single step from what it wrote.

The same modules serve both entry scripts because every step's workflow takes its input in
either shape: a samplesheet, which is what `steps.nf` passes, or a channel already carrying
one tuple per sample, which is what the previous step emits inside `main.nf`. Neither entry
script needs a module of its own.

`main.nf` renders `notebooks/celltype_report.qmd` and does not read `--qmd`; only the
`steps.nf` render steps take one.

`nextflow run .` resolves to `main.nf`, so the default answer to "run this repo" is an
analysis rather than a usage error.

### Profiles

`nextflow.config` itself carries only what every site shares — params, `workDir`, the
`--no-home`/`TMPDIR` containerOptions fix — and states no container, executor, or resource
caps of its own. Every site-specific value, including which container each step runs in,
lives in `conf/<profile>.config`. There is deliberately no `standard` profile: Nextflow only
auto-applies that when `-profile` is omitted entirely, so naming one would mean `-profile wsl`
silently dropped the defaults with it.

| Invocation | What you get |
|---|---|
| *(none)* | no container configured — steps fail rather than silently running in the wrong image |
| `-profile wsl` | local executor, the containers, how a `label 'gpu'` process reaches a card under WSL2 |
| `-profile oscer` | SLURM, scratch and OURdisk paths, the GPU queue |

Forgetting `-profile` is therefore a hard failure, not a silent login-node run — the cost
of that safety net was relying on `apptainer` config merging correctly between
`nextflow.config` and the profile, which is now avoided entirely.

Add `-stub` to check wiring without doing the work.

### How artifacts are named

Every artifact is `<sample>.<step>.<ext>` — `u2os_test.cluster_spatialdata_gpu.zarr`, not
`u2os_test.zarr`. Nothing then collides when files are staged flat: not two samples, and not
one step's input with its own output.

That is what lets a report notebook ask for a specific step —
`glob("*.annotate_celltypes.zarr")` — instead of `*.zarr` and hoping, and it is why the render
step can stage everything into one directory rather than numbered `input*/` subdirectories.

A `--group_by` run puts the column in the name — `<sample>.<column>.centroids.h5ad`, not
`<sample>.centroids.<column>.h5ad` — so one glob takes both.

Step 1a is the one exception: its three files are named for what `merscope()` looks for, not
for the sample and step. They are read by name from a directory, never staged flat, so
there is nothing for them to collide with.

Output publishes beside the code, so a result sits next to the analysis that produced it.
Each invocation gets its own directory under `results/`, named by a timestamp `run_id`:

```
<repo>/results/<run_id>/<sample>/<step>/    published step output, not committed
<repo>/results/<run_id>/<step>_samplesheet.csv    the handoff sheet, one per step that ran
<repo>/results/<run_id>/pipeline_info/    nextflow's own timeline and report
<repo>/reports/<name>/    one directory per render, committed and pushed
```

Artifacts nest by sample and step; the handoff sheets sit at the top of the run directory,
being one file per step rather than one per sample. A sheet's path column names the
published location, not the work dir, so a step re-run from a sheet reads a file that still
exists after the run that wrote it is gone.

The work dir and the image cache stay out of the repo, being large, churny, and
reproducible: under `~/merfish-mouse-testis_work/` by default, and
`/scratch/$USER/merfish-mouse-testis_work/` on `oscer`. Scratch deletes files 14 days after they
are created no matter how recently they were read, so nothing durable can live there — on
`oscer` the image cache sits on OURdisk for that reason.

The work dir is dropped when a run succeeds, so `-resume` works after a failure and not
after a success. Because runs never share an output directory, a `-stub` run cannot
overwrite a real run's results; pass `--run_id <name>` to pin the directory, which
`-resume` needs across launches.

On `oscer` the repo itself lives on OURdisk, which is permanent and large but **never
backed up**. The code and `reports/` are safe because they are pushed to GitHub, and
`results/` is reproducible from them; `data/raw/` is neither, and is worth a second copy.

## Workflow

Steps run in this order. This table is the ordering contract — there are no numeric
filename prefixes.

| Step | Script | Samplesheet | Input | Output |
|------|--------|-------------|-------|--------|
| 1 | `bin/create_spatialdata.py` | `sample, path` | MERSCOPE region directory | `<sample>.<step>.zarr` |
| 1a | `bin/prep_cellpose_vpt.py` | `sample, path, cellpose_path` | MERSCOPE region directory and a merged bespoke cellpose segmentation | `cellpose_*.{csv,parquet}` and `detected_transcripts.csv` |
| 1b | `bin/create_spatialdata_cellpose.py` | `sample, path, vpt_path` | MERSCOPE region directory and VPT cellpose output, from a VPT run or from step 1a | `<sample>.<step>.zarr` |
| 2 | `bin/cluster_spatialdata_gpu.py` | `sample, path` | zarr from step 1 or 1b | `<sample>.<step>.zarr` |
| 3 | `bin/annotate_celltypes.py` | `sample, path` | zarr from step 2 | `<sample>.<step>.zarr` |
| 4 | `bin/create_centroids.py` | `sample, path` | zarr from step 2 or 3 | `<sample>.centroids.h5ad` |
| 5 | `--qmd`, e.g. `notebooks/celltype_report.qmd` | `sample`, plus whatever path columns the .qmd globs | for `celltype_report.qmd`: zarr from step 3 and centroids from step 4 | `reports/<qmd>_<run_id>_<to>/`, one directory per render — `render_sample` nests one per sample inside it |

### 1. create_spatialdata

Reads a raw MERSCOPE output directory and writes it as a SpatialData Zarr store. The
sample id is written into `table.obs["sample"]`, which later steps read instead of
parsing a staged filename, and the loaded z-plane into `table.uns["z_layer"]` as
provenance — nothing reads that one yet.

The store holds the mosaic image as a multiscale pyramid, the transcripts as 3D points,
the cell boundaries as polygons, and the count matrix as the table.

Elements are named `<sample>_<region dir>_<element>`, from the sample id and the raw
region directory's own name — never from the staged copy below, so staging a sample
leaves its element names exactly where reading it directly would.

Two things about a region directory are fixed before the read, into a staged copy that
symlinks everything it is not changing so the instrument output is never touched. A
`cell_metadata.csv` written in a different row order than `cell_by_gene.csv` is
reordered, since the reader hands both straight to AnnData. And boundaries written the
pre-VPT way — a `cell_boundaries/` directory of per-FOV HDF5 files rather than
`cell_boundaries.parquet` — are converted, at z-index 0, the only plane the reader keeps.

A region with neither form of boundary is an error here. The reader only warns and loads
no polygons, leaving a store that clusters and annotates perfectly well and has nothing to
draw a tissue figure from, which is not a thing to discover in step 5.

The script is a plain CLI and runs outside Nextflow unchanged:

```bash
apptainer exec oras://ghcr.io/brent-biddy/python_cpu-sif:1.0.2 \
    python bin/create_spatialdata.py \
        --sample testis_01 \
        --path data/raw/testis_01 \
        --outdir results/testis_01/create_spatialdata
```

### 1a. prep_cellpose_vpt

Converts a cellpose segmentation that was **not** run through VPT into the three files
step 1b reads, so a store built from it comes out of the same reader as every other
sample. Only needed for the lab's own cellpose pipeline — `area_cellpose.py` per FOV,
`merge.py` to stitch — which writes numpy arrays rather than CSVs and a parquet.

It reads the merged segmentation directory:

```
labels.npy                    (7, height, width) uint32 raw labels, full mosaic pixels
<sample>-cellp-label-map.npz  raw label -> compact cell id, indexed by label - 1
<sample>-cell-by-gene.npz     counts, row 0 is background
<sample>-slicee.pkl           per-label bounding boxes, from ndi.find_objects
cell-coms.npz                 (z, y, x) centre per cell, mosaic pixels
cell-vols.npz                 voxels per cell
```

None of the index conventions are recoverable from the files themselves, so they are
documented in the script's docstring and taken from `merge.py`:

```
compact_id = label_map[raw_label - 1]    0 means the label is absent
counts     = cbg[compact_id]             row 0 is background, not a cell
centre     = cell_coms[compact_id - 1]
volume     = cell_vols[compact_id - 1]
```

**`label_map` is indexed by `label - 1`, not by the label.** `merge.py` builds it over
`ndi.find_objects` output, whose element *i* describes label *i+1*. Indexing by the label
happens to agree wherever labels are contiguous and is off by one across every gap, which
is silent — it gives a cell its neighbour's counts. The step therefore checks the join
before writing anything: each label's voxels in the raster must equal the volume
`merge.py` recorded for the cell it maps to, exactly, or the run fails.

The segmentation is 3D over seven z planes and a store holds one, so **each cell's
boundary is taken from the plane it covers most**. Fixing a single plane instead would
leave the ~3.5% of cells centred near the top or bottom of the section with counts and a
centre but no polygon. The cost is that a tissue figure mixes depths rather than showing
one optical section; which plane each cell came from is written to the metadata and lands
in the table as `obs["z_plane"]`. Every polygon is written with `ZIndex` 0 regardless,
because that is the only value `merscope()` keeps.

Coordinates are converted with the region's own `micron_to_mosaic_pixel_transform.csv`
rather than a restated pixel size. Gene names come from the region's `cell_by_gene.csv`
header, since the counts array is bare integers — that the blank codewords carry ~20x
lower counts than real genes is what confirms the order applies.

```bash
nextflow run steps.nf -profile oscer --step prep_cellpose_vpt \
    --samplesheet assets/samplesheet_cellpose_prep.csv

nextflow run steps.nf -profile oscer --step create_spatialdata_cellpose \
    --samplesheet results/<run_id>/prep_cellpose_vpt_samplesheet.csv
```

It writes a fourth file VPT itself only writes on request: the region's transcripts with a
`cell_id` column, `-1` where a transcript is in no cell. Assignment is the same raster
lookup that produced the counts — a transcript's `global_z` is a plane index, not a depth in
microns, so it selects the labels plane directly. Because a cell's boundary is one plane and
its transcripts are three-dimensional, a transcript can belong to a cell and still fall
outside the polygon drawn for it.

Every cell the segmentation found gets a row in the counts and the metadata, including one
whose mask traced to no polygon; `has_boundary` in the metadata says which. Only the parquet
can lose a cell, and a store missing a boundary is undrawable where a store missing the cell
is unrecoverable.

The handoff sheet it writes is already in step 1b's `sample, path, vpt_path` shape,
forwarding the region unchanged and naming this step's output as the VPT directory.

This is the memory-hungriest step: the label raster is about 7 GB per plane and the step
holds one plane, the counts and the traced polygons at once, measured at a peak near 20 GB.
It has no per-step override, so on `oscer` it takes the retry ladder — 32 GB on
the first attempt and 32 GB more on each of the three retries — and locally it takes the
16 GB default, which is not enough. It has never been run on `oscer`; the first run that
does is worth recording the real peak from.

### 1b. create_spatialdata_cellpose

An alternative to step 1 for a region that has been re-segmented with cellpose, writing a
store of the same shape so steps 2 onward read it unchanged. Either step produces a `.zarr`
and a handoff sheet; a run uses one or the other, not both.

`--vpt_path` takes either a real VPT output directory or step 1a's output, which is
written in the same shape for exactly that reason.

It reads two directories. A re-segmentation replaces the cells and which cell each
transcript belongs to, so the MERSCOPE region directory supplies only the mosaic images; the
VPT output directory supplies the rest:

```
<region>/Cellpose/cellpose_cell_by_gene.csv
                  cellpose_cell_metadata.csv
                  cellpose_micron_space.parquet
                  detected_transcripts.csv
```

`merscope()` reads the first three natively through its `vpt_outputs` argument, so nothing
is converted. They are named explicitly rather than by pointing the reader at the directory,
because only the boundaries carry the segmentation method in their name: these are the names
a Vizgen cellpose delivery uses, while a stock VPT run writes the two CSVs unprefixed and a
watershed run names its boundaries `watershed_micron_space.parquet`. Any of the three
missing is an error, naming both alternatives.

The transcripts are not among them: `vpt_outputs` has no slot for them, and `merscope()`
takes them from the region directory. VPT's copy is reachable only by staging it in over the
region's, which is what this step does — and a VPT directory without one is an error too,
since the store's points would carry no assignment while its cells are cellpose cells. A
stock VPT run only writes that file when asked; step 1a always does.

That the boundaries are an error rather than a warning is for the reason it is in step 1:
the reader only warns, and a store with no polygons is not a thing to discover in step 5.
The polygons themselves are already in micron space, so unlike the pre-VPT HDF5 form step 1
converts, there is nothing to build.

It fixes the same thing step 1 does, for the same reason: a `cellpose_cell_metadata.csv` in
a different row order than the counts is written out reordered, since `merscope()` hands both
straight to AnnData. Like step 1, the fix goes into a staged copy that symlinks everything it
is not changing, so neither the instrument output nor the VPT output is modified.

Which segmentation a store came from is not recoverable from its contents, so
`table.uns["segmentation"]` and `table.uns["vpt_path"]` record it. Give the two stores
different sample ids if you want to keep both — they are otherwise indistinguishable in a
report, and element names are built from the sample id.

```bash
nextflow run steps.nf -profile wsl --step create_spatialdata_cellpose \
    --samplesheet assets/samplesheet_cellpose.csv
```

### 2. cluster_spatialdata_gpu

Filters cells, normalizes, runs PCA / neighbors / UMAP, and sweeps Leiden resolutions on
the GPU with rapids-singlecell. The sweep is 0.1 to 2.0 in steps of 0.1, set in the script.

Each resolution leaves two obs columns, the resolution written to two decimals.
`leiden_res_0.10_v0` is Leiden's own labelling, kept so a cluster can be traced back;
`leiden_res_0.10_v1` renumbers those clusters `1..k` by descending cell count, so cluster
1 is always the largest. **v1 is what everything downstream means by a cluster id** —
Leiden's own numbers say nothing about size and are not comparable between two
resolutions, so a cluster named in a figure or a hand-written annotation is named by its
v1 id.

Cells are filtered on transcript count: `--min_counts` drops the low tail, 20 by default,
and `--max_counts_quantile` cuts doublets and segmentation merges off the top, disabled by
default. Both are recorded into `uns`, and the quantile is taken before the low cut so the
threshold does not move with it. Neither is a Nextflow param — a run that wants other
values runs the script directly.

There is no highly-variable-gene selection — a MERFISH panel is a few hundred curated
markers, so every gene is used. Genes are not filtered either.

This step needs a CUDA GPU. Driver passthrough differs by profile — `oscer` uses `--nv`,
`wsl` binds the WSL2 driver directory — and is set in the site config, not the module. The
process declares `label 'gpu'`; each site answers it.

```bash
nextflow run steps.nf -profile wsl --step cluster_spatialdata_gpu \
    --samplesheet results/<run_id>/create_spatialdata_samplesheet.csv
```

The samplesheet is the handoff sheet step 1 wrote.

### 3. annotate_celltypes

Spearman-correlates every cell against each cell type in a reference centroid table,
writing one `corr_<cell type>` obs column plus `cell_type_per_cell`, the type a cell
correlates with most strongly. Genes are matched case-insensitively, so a mouse panel's
`Acta2` finds a human reference's `ACTA2`; which panel genes were found is written to
`<sample>.gene_overlap.tsv`.

Correlations are standardized within each cell, because everything downstream compares
cells with each other and a cell's correlation to every type rises with how many genes it
captured. That cannot change which type is largest, so the calls are unaffected.

**This step makes no cluster-level call.** A cluster's identity is read off the composition
of these per-cell calls, which is a judgment made by a person — the reference reports where
a type would land, not that it was found.

`--reference` selects the table; it defaults to
`assets/reference/shami_human_testis_centroids.csv.gz`, the only one there.

```bash
nextflow run steps.nf -profile wsl --step annotate_celltypes \
    --samplesheet results/<run_id>/cluster_spatialdata_gpu_samplesheet.csv
```

### 4. create_centroids

Builds one row per cluster from the clustered zarr, at every resolution in the sweep, so
later steps and reports never open the counts matrix. Writes a small h5ad holding summed
CP10K in `X` and summed raw counts in `layers["counts"]`, with `n_cells` per row.

Everything is a **sum, not a mean**, so clusters pool: any union's profile is the row-wise
sum of its members, and `n_cells` sums with it. The reference centroids in
`assets/reference` are `ln(mean + 1)`, so the comparable value here is `log1p(X / n_cells)`.

Its own step rather than part of step 2: Nextflow hashes the task script, so folding it in
would make a change to the centroid recipe re-run the GPU Leiden sweep.

`--group_by <column>` sums over one named obs column instead of the sweep — cell type,
once a later step has written it. Those runs are named for the column
(`<sample>.<column>.centroids.h5ad`, and `create_centroids_<column>_samplesheet.csv`), so
they publish beside a sweep run rather than displacing it and the step can be re-run for
each grouping you want.

The handoff sheet is `sample,path,centroid_path`: it forwards the zarr this step read
alongside the centroids it wrote, because a report wants both and only this step knows
which zarr the centroids came from.

A grouping that is a union of v1 clusters needs no run at all — sums are additive, so add
the rows.

```bash
nextflow run steps.nf -profile wsl --step create_centroids \
    --samplesheet results/<run_id>/cluster_spatialdata_gpu_samplesheet.csv

nextflow run steps.nf -profile wsl --step create_centroids --group_by cell_type \
    --samplesheet results/<run_id>/cluster_spatialdata_gpu_samplesheet.csv
```

### 5. render_cohort and render_sample

Renders the .qmd named by `--qmd`. Terminal steps — nothing consumes them.

Two steps, one workflow. The step name says how rows are grouped — `render_cohort` renders
every row in one pass, `render_sample` renders one per row — and `--to` says which of the
formats the .qmd declares comes out of it. It defaults to `pptx`, and one render is one
format, so both packagings means running twice.

They are separate because the right grouping differs by format. A markdown document is
navigated by file, so fifteen samples in one page is an enormous scroll where a directory of
fifteen pages is browsable; a deck is read linearly in a meeting, where one file with fifteen
sections works and fifteen files is miserable.

```bash
--step render_cohort --to pptx    # one deck, every sample
--step render_sample --to gfm     # one page per sample
```

The notebook needs no conditional for this: it globs what was staged, so a cohort render finds
every sample and loops over them while a per-sample render finds one and loops once.

The step knows nothing about what it is rendering. It requires only a `sample` column and
stages every other column's path; which of them a notebook wants, and what it does with
them, is the notebook's own business. **A new report is therefore a new
notebook and no Nextflow at all** — copy an existing `.qmd`, edit it, and render it.

Staging is the input contract: the workflow drops every samplesheet path flat beside the
notebook and the notebook globs the step it wants — `*.annotate_celltypes.zarr`, not
`*.zarr` — taking each sample's id from inside its object rather than from the staged
filename. Flat works because `<sample>.<step>.<ext>` already makes every file distinct, so
a step that names no columns needs no scheme to keep them apart. Adding a sample needs no
edit to the notebook.

`notebooks/celltype_report.qmd` is the report this pipeline has today: a cohort deck and a
GitHub-readable document with a section per sample — its QC, the clustering, the per-cell
calls, what they compose to per cluster, and the calls on tissue. The rest of this section
describes that notebook rather than the step.

It reads step 4's handoff sheet, which forwards both the zarr it consumed and the
centroids it wrote. The centroids are what the cluster similarity figure correlates, so
the notebook never opens a counts matrix.

Two knobs live at the top of the notebook, both judgments rather than computations. The
resolution the whole document reads, since the sweep writes twenty and one has to be
picked; and the share of a cluster's cells one type must take for the cluster to be
settled on it. A cluster below that is called Ambiguous rather than named on a split
vote, and gets a slide of its own laying out the evidence.

Clusters carry a second, render-only numbering: `v2` is position down the dendrogram, so
neighbouring numbers are transcriptionally adjacent and the heatmaps read as a block
diagonal. Figures are labelled `v1 (v2)`. **`v1` is what a call is written against** — it
is fixed in the object, while `v2` moves whenever the clustering or the gene set changes.

These are the only steps that publish into the repo. They write `reports/`, not
`results/<run_id>/`, because their product is the one that is source-controlled: the
markdown and its figures are committed, being what GitHub renders, and the deck is
gitignored, being the same content in a format git cannot diff.

Each render gets its own directory, `reports/<qmd>_<run_id>_<to>/`, holding the
document, the figure directory and the deck together. Quarto writes all three into
`--output-dir` and the document's links are relative to it, so the directory moves as a unit
and renders on GitHub wherever it sits. The .qmd leads the name so a listing groups a
report type together; the run id names the render invocation, not the results tree it read —
those are different runs, since a render takes a samplesheet pointing at some earlier run's
published paths, and what a render was built from is in that samplesheet. The format is part
of it because the two packagings of one notebook are different artifacts produced by separate
runs, and sharing a name put the per-sample pages inside the cohort render's own directory;
it is always appended, one format per render -- `--to` defaults to `pptx`, so rendering both
packagings means running twice. A `render_sample` run nests one directory per sample inside
that. `--run_id` names a render when it deserves a name of its own:

```bash
# One deck over the cohort.
nextflow run steps.nf -profile wsl --step render_cohort --to pptx --run_id cellpose_cmp \
    --qmd notebooks/celltype_report.qmd \
    --samplesheet results/<earlier_run>/create_centroids_samplesheet.csv
# -> reports/celltype_report_cellpose_cmp_pptx/celltype_report.pptx

# The same samples as one browsable page each.
nextflow run steps.nf -profile wsl --step render_sample --to gfm --run_id cellpose_cmp \
    --qmd notebooks/celltype_report.qmd \
    --samplesheet results/<earlier_run>/create_centroids_samplesheet.csv
# -> reports/celltype_report_cellpose_cmp_gfm/<sample>/celltype_report.md
```

A render whose prose is specific to one comparison is a copied `.qmd` rather than a note
injected from outside, because commentary about a figure belongs next to the figure. The
cost is that a fix to a shared figure has to land in each copy.

One collision is left on purpose: the same format in both modes — `render_cohort --to gfm`
alongside `render_sample --to gfm` — still shares a directory, since the grouping is not in
the name. The pairing above avoids it; making it impossible would mean naming the grouping
too, for a case nothing needs yet.

Apart from that one case, no two renders share a directory, so nothing overwrites anything
— which also means `reports/` accumulates. Prune the ones not worth keeping before
committing.

`reports/README.md` indexes the ones that were kept, and is what GitHub shows when anyone
browses the directory. It is written by hand: a render cannot say what it was for, and a
pipeline that wrote the file would overwrite the part worth reading. Add a row when you
commit a report.

Figure resolution is set per format in the notebook's frontmatter — 144 dpi for the
committed markdown, read in a ~900px column, and 200 for the projected deck. It cannot go
back into `rcParams`, which runs after quarto and would override both.

Two segmentations of one sample are compared by giving them distinct sample ids —
`testis_01_vizgen` and `testis_01_cellpose` — in one samplesheet. Both readers write the
id into `table.obs["sample"]`, and the notebook takes each section's id from inside the
object, so the two land as neighbouring sections of a single document with no edit to the
notebook and no collision in `results/`.

`error: true` is not set, so an uncaught failure fails the render and its exit code. The
sample loop catches its own, though: a store that will not open, one missing an obs column
an earlier step should have written, or one slide that raises all draw into the document
and the other samples still render. So the exit code answers "did setup work", not "is the
deck complete" — inspect it either way: `unzip -q reports/<dir>/celltype_report.pptx` and
check the slide count and titles.

## Layout

```
main.nf          a chained analysis; `nextflow run .` resolves here
steps.nf         step dispatch, for working one step at a time
nextflow.config  params, the defaults that always apply, and which site profiles exist
conf/            one file per site: wsl and oscer state their difference
bin/             all executable code
modules/         one file per step: its process and its workflow
modules/samplesheet.nf   the one shared helper: how a step reads its input
notebooks/       report notebooks
assets/          sample sheets, and the pptx template and lua filter a render needs
assets/reference/  cell type centroids to annotate against; see each file's header
data/raw/        raw instrument output (not committed)
results/<run_id>/<sample>/<step>/    published step output (not committed)
results/<run_id>/<step>_samplesheet.csv    the handoff sheet each step writes
reports/README.md  hand-written index of the renders worth keeping
reports/<qmd>_<run_id>_<to>/    one directory per render: markdown and figures
                 committed, the deck not; render_sample nests one dir per sample
```

`reports/` is the only output that is committed. Everything else under the repo is
gitignored and reproducible from `bin/` + `assets/` + `data/raw/`.
