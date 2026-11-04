# Notebook design

How the analysis notebook is wired. `PLAN.md` lists its steps in order and what each one
decides; this file is how they are coded: the calls and column names for each step, and what
cuts across them: the one object, the checkpoints, the reference cache, and CPU or GPU. The
recipe is the earlier pipeline's, carried over, on scanpy (`sc`).

The analysis is one self-contained Quarto notebook, run once per section. It calls no `bin/`
scripts: everything from building the reference to the cluster calls is code in the notebook,
built up a section at a time in the order of `PLAN.md`'s steps. A step is named by phase and
number, B1 or D5. Its input is the store and table Phase A (pipeline steps 1a and 1b) already
built for the section.

It is worked out end to end on one section, b2r0, before it runs on the other 15: every phase
up to I is settled there first, and Phase J follows once every section is through I. Anything
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
| Reference | D0 | The reference and its level; no sample, so every section shares it |
| 1, embedded | C6 | The section, the filter, normalization, PCA and neighbour settings, and the backend |
| 2, annotated | D8 | Checkpoint 1's key, the resolution, the reference and the purity threshold |

Changing anything in Phase C reruns C1 to C7 and invalidates checkpoint 2 with it, since
checkpoint 2's key includes checkpoint 1's. Changing anything in Phase D reruns only D1 to D8
from checkpoint 1, so annotation is iterated on without reclustering, and a run with each
reference leaves its own checkpoint 2 side by side. Checkpoint 2 is the whole object, so the
later phases read whichever columns they need from it.

A checkpoint after C2 could be added so that changing PCA, neighbours or Leiden does not redo
filtering and normalization. Not needed yet; add it if clustering b2r0 turns out slow.

## One clustering function, two levels (deferred)

Steps C2 to C7 are one function, `cluster_cells(adata)`, and the same function reclusters a
single cell type's cells in D9. Nothing is reused between levels: the subset starts again from
its raw counts layer and runs the same normalize, scale, PCA, neighbours, UMAP and Leiden calls
with the same parameters. Only the filter (C1) is skipped, since those cells already passed it.

- Normalization means the subset's own median total, so log-normalized values are on a
  different scale at each level. Anything compared across levels starts from counts.
- A subset too small for 30 PCs or k = 15 raises an error rather than quietly changing the
  recipe.
- Four resolutions are written at each level so annotation can move to another without
  reclustering.

## Step code and column names

Keyed by the step numbers in `PLAN.md`. The per-cell method in Phase D is the earlier
pipeline's, carried over.

**Phase C**

| Step | Code |
|---|---|
| B1 | `adata = anndata.read_h5ad(TABLE)`, where `TABLE` is the one `*.create_spatialdata_cellpose.h5ad` staged beside the notebook; `adata.layers["counts"] = adata.X.copy()`; then `adata.obs["transcript_count"] = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()` and `adata.obs["gene_count"] = np.asarray((adata.layers["counts"] > 0).sum(axis=1)).ravel()`. The table has neither column, so B1 derives both. They cover the 491 real genes only: blank codewords are in `obsm["blank"]`, not in `X`, and are left out. The sample id is read from `adata.obs["sample"]`, not from the file name. Phase B's plots and C1's high cut read `transcript_count`. B1 also writes `adata.obs["table_row"] = np.arange(adata.n_obs)` and keeps `N_TABLE_CELLS = adata.n_obs`, because C1 drops cells and an external per-cell file is indexed by the full table |
| C1 | `sc.pp.filter_cells(adata, min_counts=20)`. The high cut's threshold is `adata.obs["transcript_count"].quantile(q)`, taken before the low cut, then `sc.pp.filter_cells(adata, max_counts=threshold)` |
| C2 | `sc.pp.normalize_total(adata)`, `sc.pp.log1p(adata)`; then `adata.layers["lognorm"] = adata.X.copy()` |
| C3 | `sc.pp.scale(adata, zero_center=False, max_value=10)` |
| C4 | `sc.pp.pca(adata, n_comps=30, random_state=0)` |
| C5 | `sc.pp.neighbors(adata, n_neighbors=15, metric="cosine", random_state=0)` |
| C6 | `sc.tl.umap(adata, random_state=0)` |
| C7 | *Not run; the notebook does not cluster.* For each resolution `r`: `sc.tl.leiden(adata, resolution=r, key_added=f"leiden_{r}", flavor="igraph", n_iterations=2, random_state=0)`, then our own `relabel_by_size()` renumbers the clusters 1..k by size |

**Phase D**

| Step | Code |
|---|---|
| D0 | Built per reference as described under the reference cache below; the chosen one is stored in `adata.uns["reference"]` |
| D1 | *Deferred.* `cluster_key = "leiden_1.0"` |
| D2 | `sc.tl.rank_genes_groups(adata, groupby=cluster_key, layer="lognorm", method="wilcoxon")`, read back with `sc.get.rank_genes_groups_df` |
| D3 | `scipy.stats.rankdata(adata[:, shared].layers["counts"], axis=1)` for the cells, the same for each reference type, then `1 - scipy.spatial.distance.cdist(cell_ranks, type_ranks, metric="correlation")`. `shared` is matched case-insensitively |
| D4 | Subtract each cell's mean correlation and divide by its standard deviation; write `adata.obs[f"corr_{type}"]` for every type |
| D5 | `idxmax(axis=1)` over the standardized correlations, written to `adata.obs["cell_type"]`; then `adata.obs["germ_somatic"]` is `"germ"` where `cell_type` is in `GERM = ["SPG", "Scytes", "STids", "Elongating"]` and `"somatic"` otherwise. This is the first pass's annotation |
| D6 | *Deferred.* `pd.crosstab(adata.obs[cluster_key], adata.obs["cell_type"], normalize="index")`, kept in `adata.uns`; `germ_share = composition[GERM].sum(axis=1)`, where `GERM = ["SPG", "Scytes", "STids", "Elongating"]` |
| D7 | *Deferred.* `cluster_type` is `"Germ"` where `germ_share` is at least `PURITY` (0.90), even when one germ type alone is; otherwise the top type where its share is at least `PURITY`, and `"Ambiguous"` where neither is, written to `adata.obs["cluster_type"]`; how it was made to `adata.obs["cluster_call_source"]`, one of `auto`, `manual` or `unresolved`. `cell_type` is left as the per-cell call, with all four germ types |
| D8 | *Deferred.* A manual call is written to `cluster_type` with source `manual`. A cluster left unresolved keeps `Ambiguous` with source `unresolved`, and its cells get `adata.obs["excluded"] = True` and `adata.obs["excluded_at"] = "major"` |
| D9 | `cluster_cells` on one major type's cells, then D3 to D5 on the result (D6 to D8 too, if turned on). For germ cells that is `adata[adata.obs["germ_somatic"] == "germ"]`, restarted from `layers["counts"]`. The output is `adata.obs["fine_cell_type"]`. A type with no fine reference keeps its major label there, with `adata.obs["fine_source"] = "major"`. An unresolved subcluster's cells get `excluded_at = "fine"` |

**Phase E**

| Step | Code |
|---|---|
| E1 | `cell_map = np.load(TUBULE_MAP, allow_pickle=True)["cell_map"]`, where `TUBULE_MAP` is the one `*cell-to-tubule-map-goodcells.npz` staged beside the notebook; `assert len(cell_map) == N_TABLE_CELLS`; `adata.obs["tubule_id"]` is `cell_map[adata.obs["table_row"]]` as a string, `""` where the cell has no tubule; `adata.obs["in_tubule"] = adata.obs["tubule_id"] != ""`. Tubule ids are strings because some sections' ids are not integers |
| E2 | `tubule_meta = pd.read_csv(TUBULE_META, index_col=0)`, the one `*tubule-metadata*.csv` staged (a `-staged` file for some sections), its index as strings. `tubules` is a DataFrame of `tubule_meta` plus `n_cells` and `germ_share` per tubule from `adata.obs.groupby("tubule_id")`, with columns `good_cells`, `good_circle`, `good_cells_bin`, `area`, `eccen`, `maj_axis_length`. Areas and lengths are in mosaic pixels |
| E3 | `plot_tubules_on_tissue` and `plot_tubule_qc`; the share inside by cell type is `adata.obs.groupby("cell_type", observed=True)["in_tubule"].mean()` |

**Phase G**

| Step | Code |
|---|---|
| G1 | `STAGED = "stage" in tubules.columns`, true when the metadata file is a `-staged` one. `staged` is the rows of `tubules` whose `stage` is in `STAGE_GROUPS = ["1-3", "4-6", "7-8", "9-12"]` and whose `good_cells_bin` is true, with `stage` as an ordered categorical. Finer values such as `3-4` or `7` are left out and counted |

**Phase H**

| Step | Code |
|---|---|
| H1 | `composition = pd.crosstab(inside["tubule_id"], inside["cell_type"].astype(str), normalize="index").reindex(columns=CELL_TYPES, fill_value=0).loc[staged.index]`; the mean per stage group is `composition.groupby(staged["stage"], observed=True).mean()`; `scipy.stats.kruskal` per type across the four groups. `inside` is the cells with `in_tubule` |
| H2 | `counts` is the same crosstab without `normalize`; `area_um2 = staged["area"] * UM_PER_PX ** 2` (`UM_PER_PX = 0.108`); `density = counts.div(area_um2, axis=0) * 1e4`, cells per 10,000 µm², with `germ` and `all cells` columns added |

Phases E to J read `fine_cell_type` and drop cells where `excluded` is set. In Phase E, per-cell
results (`tubule_id`, `in_tubule`) go in the table. The tubule polygons are not loaded yet; when
Phase F needs them they are Alex's `global_path` arrays, saved beside the table.

## The reference cache

D0 is the first step any run reaches that needs the reference. Its checkpoint is written at
`results/<notebook>/reference/` rather than under a sample, so it is made once outside the
notebook and all sixteen sections load a table of a few hundred genes by eleven types. The full GEO dataset
is read in chunks and dropped once the table is saved, so it never stays in memory. The
reference is also stored in the object's `uns`, so a checkpoint records which reference its
calls were made against.

How each one is made, outside the notebook:

- *Green, 11 major types:* read the GEO counts in chunks, `pd.read_csv(dge, sep="\t",
  index_col=0, chunksize=...)`, summing every cell's total over all genes and keeping only the
  panel genes' rows. Read the per-cell types with `pd.read_csv(attributes, sep="\t",
  comment="#", header=None)`. Then TP10K = panel counts / cell total x 10,000, mean per type
  with `groupby(...).mean()`, and `np.log1p`.
- *Green, GC1-12:* `pd.read_excel(xlsx, sheet_name="Average")`, the panel genes' rows.
- *Shami:* `pd.read_csv` of the table already in the repo.
- Panel genes renamed since the reference was made are mapped to the names it uses, with our
  own `match_panel_genes()`, rather than dropped.

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
