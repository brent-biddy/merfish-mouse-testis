# Notebook design

How the analysis notebook is wired. `PLAN.md` lists its steps in order, with the functions each
one calls; this file is what cuts across those steps: the one object, the checkpoints, the
reference cache, and CPU or GPU.

The analysis is one self-contained Quarto notebook, run once per section. It calls no `bin/`
scripts: everything from building the reference to the cluster calls is code in the notebook,
built up a section at a time in the order of `PLAN.md`'s steps. A step is named by phase and
number, B3 or C5. Its input is the store and table the segmentation steps (1a and 1b) already
wrote for the section.

It is worked out end to end on one section, b2r0, before it runs on the other 15: every phase
up to H is settled there first, and Phase I follows once every section is through H. Anything
learned on a later section that changes a step goes back to b2r0 too, so all 16 are analysed
the same way. `NOTES.md` tracks where the runs stand.

## One object

**One object passes through the whole notebook.** Step B1 loads the section's table as a
single AnnData, and every later section adds to that same object rather than making a new one:
raw counts kept in `layers["counts"]` and log-normalized values in `layers["lognorm"]`, since
normalization and scaling overwrite `X`, embeddings in
`obsm`, clusters, scores and calls in `obs`, markers and parameters in `uns`. A checkpoint is
that object written to disk, and loading one puts the notebook back where it was. The
SpatialData store's images and boundaries are not loaded into it; tissue figures read the
shapes they need from the store when they draw.

Memory clean-up is not planned. Adding to one object in place leaves no stale copies to free,
and the object is small: about 147,000 cells by 491 genes is under 300 MB dense, plus the
neighbour graph. Two habits keep it that way: no `adata.copy()` per section, and anything big
that is only needed briefly, such as a full reference dataset, is deleted once used.

## Checkpoints

Each section caches its result as an on-disk checkpoint: a file (`.h5ad`, `.parquet`) named
by a key built from that section's parameters and the key of the section before it. A section
whose file exists loads it instead of recomputing, so editing a later section never reruns an
earlier one, and changing a parameter reruns only that section and the ones after it. A code
change the parameters do not capture needs the checkpoint deleted, or the section's code
hashed into its key. The checkpoints are ordinary files, readable outside the notebook.

The checkpoints are the notebook's intermediate results, written into `results/` like the
pipeline's and gitignored with it:

```
results/<notebook>/<sample>/<section>/<sample>.<section>.<key>.<ext>
```

The layout is the pipeline's `results/<run_id>/<sample>/<step>/` with the notebook's name in
place of the run id. A timestamped run id would start every render in an empty directory and
find nothing to load; a fixed one lets each render reuse what the last one wrote. The key in
the filename keeps results from different parameters side by side rather than overwriting
each other.

## Where the checkpoints fall

| Checkpoint | After | Key holds |
|---|---|---|
| Reference | C0 | The reference and its level; no sample, so every section shares it |
| 1, clustered | B8 | The section, the filter, normalization, PCA, neighbour and Leiden settings, and the backend |
| 2, annotated | C8 | Checkpoint 1's key, the resolution, the reference and the purity threshold |

Changing anything in Phase B reruns B2 to B8 and invalidates checkpoint 2 with it, since
checkpoint 2's key includes checkpoint 1's. Changing anything in Phase C reruns only C1 to C8
from checkpoint 1, so annotation is iterated on without reclustering, and a run with each
reference leaves its own checkpoint 2 side by side. Checkpoint 2 is the whole object, so the
later phases read whichever columns they need from it.

A checkpoint after B3 could be added so that changing PCA, neighbours or Leiden does not redo
filtering and normalization. Not needed yet; add it if clustering b2r0 turns out slow.

## One clustering function, two levels

Steps B3 to B8 are one function, `cluster_cells(adata)`, and the same function reclusters a
single cell type's cells in C9. Nothing is reused between levels: the subset starts again from
its raw counts layer and runs the same normalize, scale, PCA, neighbours, UMAP and Leiden calls
with the same parameters. Only the filter (B2) is skipped, since those cells already passed it.

- Normalization means the subset's own median total, so log-normalized values are on a
  different scale at each level. Anything compared across levels starts from counts.
- A subset too small for 30 PCs or k = 15 raises an error rather than quietly changing the
  recipe.
- Four resolutions are written at each level so annotation can move to another without
  reclustering.

## The reference cache

C0 is the first step any run reaches that needs the reference. Its checkpoint is written at
`results/<notebook>/reference/` rather than under a sample, so the first section run builds it
and the other fifteen load a table of a few hundred genes by eleven types. The full GEO dataset
is read in chunks and dropped once the table is saved, so it never stays in memory. The
reference is also stored in the object's `uns`, so a checkpoint records which reference its
calls were made against.

The GEO files, as checked on 2026-11-03:

| File | Size | Holds | Used for |
|---|---|---|---|
| `GSE112393_MergedAdultMouseST25_DGE.txt.gz` | 81 MB | Raw UMI counts, 37,241 genes by 34,633 cells | The 11 major types |
| `GSE112393_MergedAdultMouseST25_PerCellAttributes.txt.gz` | 1.5 MB | Per cell: dataset, `CellType` as one of the 11 major types, QC metrics | The 11 major types |
| `GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx` | 11 MB | The authors' centroids for GC1-12 over 24,475 genes | GC1-12, taken as published |

Computed the paper's way, the major types match the collaborator's somatic columns exactly, and
GC1-12 match its `Germ_1` to `Germ_12`, so the collaborator's file and the public data agree
wherever both exist. GEO has no per-cell labels below the 11 major types, so the SPG and
Sertoli subtypes cannot be rebuilt. Seven panel genes are missing from the counts under the
panel's names, and nine from the germ-cell spreadsheet, being newer symbols for renamed genes
(`Prxl2a`, `Plpp3`, `Mesd`, `Meioc`, ...).

## CPU or GPU

**The notebook starts on CPU (scanpy).** b2r0 has about 147,000 cells after filtering, over
491 genes, which CPU clusters in minutes, and a CPU job runs on any OSCER node or locally with
little or no wait.

GPU (rapids-singlecell, as the earlier pipeline used) is the option if clustering becomes the
slow part. It clusters faster, but OSCER's GPU queue is the wait: the earlier pipeline's GPU
jobs on `sooner_gpu_test` were scheduled hours apart, the last a day out. The `wsl` machine's
own GPU has no queue. rapids-singlecell mirrors scanpy's API (`rsc.pp`, `rsc.tl` for `sc.pp`,
`sc.tl`), so moving over is a switch in the notebook rather than a rewrite. The two Leiden implementations do not give identical clusters, so the backend is
part of checkpoint 1's key.
