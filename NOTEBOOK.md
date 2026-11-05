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
through Phase H is settled there first, and the further analyses follow once every section is through it. Anything
learned on a later section that changes a step goes back to b2r0 too, so all 16 are analysed
the same way. The tracking issues record where the runs stand.

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
single cell type's cells later. Nothing is reused between levels: the subset starts again from
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

**Phases B and C**

| Step | Code |
|---|---|
| B1 | `adata = ad.read_h5ad(TABLE)`, where `TABLE` is the one `*.create_spatialdata_cellpose.h5ad` staged beside the notebook; `adata.layers["counts"] = adata.X.copy()`. The sample id is read from `adata.obs["sample"]`, not from the file name, and the table must hold one sample. B1 also writes `adata.obs["table_row"] = np.arange(adata.n_obs)` and keeps `N_TABLE_CELLS = adata.n_obs`, because C1 drops cells and an external per-cell file is indexed by the full table. The notebook does not yet derive `transcript_count` or `gene_count`; the QC plots that would read them, and C1's high cut, are not built. They would cover the 491 real genes only: blank codewords are in `obsm["blank"]`, not in `X` |
| B2 | `sdata = sd.read_zarr(STORE)`; the first pyramid level whose longest side is at most 6,000 pixels, `.sel(c="DAPI")`; `um_to_px` from `get_transformation(sdata.shapes[cells_name], to_coordinate_system="global")`, inverted to `px_to_um`, which also puts the image's extent in micrometres. Kept as `TISSUE_IMAGE = (dapi, extent, (low, high))` |
| C1 | `sc.pp.filter_cells(adata, min_counts=20)`. The high cut's threshold is `adata.obs["transcript_count"].quantile(q)`, taken before the low cut, then `sc.pp.filter_cells(adata, max_counts=threshold)` |
| C2 | `sc.pp.normalize_total(adata)`, `sc.pp.log1p(adata)`; then `adata.layers["lognorm"] = adata.X.copy()` |
| C3 | `sc.pp.scale(adata, zero_center=False, max_value=10)` |
| C4 | `sc.pp.pca(adata, n_comps=30, random_state=0)` |
| C5 | `sc.pp.neighbors(adata, n_neighbors=15, metric="cosine", random_state=0)` |
| C6 | `sc.tl.umap(adata, random_state=0)` |

**Phase D**

| Step | Code |
|---|---|
| D0 | Built per reference as described under the reference cache below; the chosen one is stored in `adata.uns["reference"]` |
| D3 | `scipy.stats.rankdata(adata[:, shared].layers["counts"], axis=1)` for the cells, the same for each reference type, then `1 - scipy.spatial.distance.cdist(cell_ranks, type_ranks, metric="correlation")`. `shared` is matched case-insensitively |
| D4 | Subtract each cell's mean correlation and divide by its standard deviation; write `adata.obs[f"corr_{type}"]` for every type |
| D5 | `idxmax(axis=1)` over the standardized correlations, written to `adata.obs["cell_type"]`; then `adata.obs["germ_somatic"]` is `"germ"` where `cell_type` is in `GERM = ["SPG", "Scytes", "STids", "Elongating"]` and `"somatic"` otherwise. This is the first pass's annotation |
| D6 | The germ cells alone, D3 to D5 again against `green_gc.csv` (`GC_TYPES`, GC1 to GC12): `is_germ = (adata.obs["germ_somatic"] == "germ").to_numpy()`, `spearman_to_reference` on `layers["counts"][is_germ]` over the shared genes, standardized within each cell, `argmax` to `adata.obs["germ_type"]` (a categorical, NaN for a cell that is not germ), and `adata.obs[f"corr_{gc}"]` for each. Nothing is reclustered. `gc_calls` is `germ_type` as a string, `"not germ"` elsewhere |

**Phase E**

| Step | Code |
|---|---|
| E1 | `tubules_um`: a dict of tubule id (string) to a shapely polygon in micrometres, `affine_transform(Polygon(np.load(path)["global_path"][:, ::-1]), to_um)`, where `to_um` is the inverted pixel transform from B2 and no shift is applied; `tubule_meta = pd.read_csv(TUBULE_META, index_col=0)`, its index as strings, the one `*tubule-metadata*.csv` staged (a `-staged` file for some sections). `flags` gives each tubule one of `TUBULE_FLAGS`: good circle, good cells but not circular, not good cells, or no metadata |
| E2 | `plot_tubule_overview`: the outlines over DAPI by flag, and over the per-cell type and the germ-cell type |
| E3 | `cell_map = np.load(TUBULE_MAP, allow_pickle=True)["cell_map"]`, where `TUBULE_MAP` is the one `*cell-to-tubule-map-goodcells.npz` staged beside the notebook; `assert len(cell_map) == N_TABLE_CELLS`; `adata.obs["tubule_id"]` is `cell_map[adata.obs["table_row"]]` as a string, `""` where the cell has no tubule; `in_tubule = adata.obs["tubule_id"] != ""`. Tubule ids are strings because some sections' ids are not integers |

**Phase F**

| Step | Code |
|---|---|
| F1 | `stage_table = pd.read_csv(TUBULE_ORDER)`, filtered to `samp == section`, with `tube_id` as a string; `stage_of = stage_table.set_index("tube_id")["prop_stage"].astype(int)`, groups 1 to 4. Also keeps `ang` and `order`. `tubule_order` is `stage_table.sort_values("order")["tube_id"]`, which is the x axis of F5. A tubule not in `stage_of` has no stage and is left out of any plot ordered by it. `ang` is circular and is not correlated with anything directly |
| F2 | `plot_stage_outlines`: the staged tubules thick in `STAGE_COLORS[stage]` and the others thin and grey, on DAPI and over the germ-cell type calls |
| F3 | `stage_calls = adata.obs["tubule_id"].map(stage_label)`: `stage N`, `unstaged tubule` or `outside tubules`; drawn with `plot_calls`, the grey groups first |
| F1 (phase) | The continuous stage, `phase`. `ring_start` is taken from the whole order file before it is filtered to the section: `ring_start = order_all.loc[order_all["order"] == 0, "ang"].iat[0]`, the angle where group 1 begins on the pooled ring. Rank rises as the angle falls, so `phase = (ring_start - stage_table["ang"]) % (2 * np.pi)`, indexed by `tube_id` like `stage_of`, running 0 to 2π. Check that `phase` rises with `order` within each section, apart from ties. `cuts` are the four group boundaries on that scale: `edges = phase.groupby(stage_of).agg(["min", "max"])`, and each of the three cuts is the midpoint between one group's `max` and the next group's `min`; the first group starts at 0 and the last ends at 2π. The rank correlation is asserted above 0.95 (it is 1.0 on all sixteen sections) |
| F2 (phase) | A third panel on `plot_stage_outlines`: the same window, DAPI beneath, the staged outlines in `plt.cm.twilight(phase / (2 * np.pi))` and the unstaged ones thin and grey as in the other two. The function takes `phase_colors` (a dict of tubule id to RGBA) beside `stage_colors`. The colourbar is `ScalarMappable(Normalize(0, 2 * np.pi), PHASE_CMAP)` with long ticks at `[0, *cuts, 2 * np.pi]` and the group names as minor ticks at the middle of each span, in `add_phase_colorbar`. The figure is `figsize=(27, 8.6)` for three panels |
| F3 (phase) | A second tissue and UMAP figure beside F3, with each cell given its tubule's phase: `cell_phase = adata.obs["tubule_id"].map(phase)`, NaN for a cell outside every tubule and for one in an unstaged tubule. The NaN cells are drawn first in grey and the rest on top with `c=cell_phase, cmap="twilight", vmin=0, vmax=2 * np.pi`, as F3 does for the groups, and the same colourbar with the cuts. A cell takes its tubule's value, so this is not independent of the calls: the stage was computed from gene profiles of these same cells |
| F1 (manual) | For the sections that have one: `STAGED = "stage" in tubules.columns`, true when the metadata file is a `-staged` one. `staged` is the rows of `tubules` whose `stage` is in `STAGE_GROUPS = ["1-3", "4-6", "7-8", "9-12"]` and whose `good_cells_bin` is true, with `stage` as an ordered categorical. Finer values such as `3-4` or `7` are left out and counted |

**Phase G**

| Step | Code |
|---|---|
| G1 | `radial_position(xy, polygon)` returns `(rad, dist_to_edge, inside)` for an `(n, 2)` array of micrometre points, after his `create_tubule_transcript_polar_info`. `c = polygon.centroid`; `angle = arctan2(y - c.y, x - c.x)`; `rad = hypot(x - c.x, y - c.y)`; `inside = shapely.contains_xy(polygon, x, y)`, so a point on the boundary is outside. For the points inside, a ray from the point to `point + reach * (cos(angle), sin(angle))`, with `reach` the length of the polygon's bounding-box diagonal, as `shapely.linestrings`; `hits = shapely.intersection(polygon, rays)`. `dist_to_edge` is `shapely.length(hits)` where the hit is a Point or a LineString, and otherwise the largest distance from the point to any coordinate of the hit's line parts. It is 0 for a point not inside. It is applied per tubule to the cells with that `tubule_id`, on `tissue_xy`, and `adata.obs["radial_pos"] = rad / (rad + dist_to_edge)`, NaN outside tubules. Equal to his `scaled_rad` exactly on b2r0; about 30 s |
| G2 | `bin_edges = np.sqrt(np.arange(N_BINS + 1) / N_BINS)`, `N_BINS = 50`. `tubule_order` is the staged tubules sorted by `phase`; `tubule_area` is `tubules_um[t].area`. `radial_bin = np.clip(np.searchsorted(bin_edges, radial_pos, side="left") - 1, 0, N_BINS - 1)`, so a bin is open at its low edge. `counts` has shape `(types, tubules, bins)` for `DENSITY_TYPES` (SPG, Scytes, STids, Elongating, Sertoli, Leydig), from `np.add.at` on `cell_type`. `density = counts / (tubule_area[None, :, None] / N_BINS) * 1e4`, cells per 10,000 µm², not rescaled per tubule. `smoothed` is `gaussian_filter(panel, sigma=RADIAL_SIGMA, mode=("wrap", "reflect"))` for each type, `RADIAL_SIGMA = (2, 1.5)`. `plot_density_heatmaps` draws both rows with `imshow(arr.T, origin="lower")`, magma, each panel from 0 to its own 99th percentile, and the stage groups above |
| G3 | `farthest_crossing`, a numba loop over a polygon's edges: for each point, the largest `t` in `(0, reach]` at which its ray `(dx, dy)` crosses an edge `a -> b` (`t = cross(a - p, e) / cross(d, e)`, `u = cross(a - p, d) / cross(d, e)`, with `0 <= u < 1`). That is the length of the one piece the ray cuts from the outline, or the farthest coordinate where it cuts several, so it equals G1's `dist_to_edge` to 3e-12 µm and runs 12 million transcripts in about 25 s. `radial_position_fast(xy, polygon)` wraps it. The transcripts are `sdata.points[...][["x", "y", "gene"]].compute()`, `gene` mapped to `adata.var_names` position, which drops the blanks. `contained[name]` is the transcripts in the outline's bounding box where `shapely.intersects_xy` is true. His assignment: `owner = np.full(n, -1)`, then for each tubule in order of `int(id)`: `at = contained[name]`, `claimed = owner[at]`; a tubule with `good_cells_bin` sets `owner[at] = k`; one in `tubule_meta` and not good sets it where `claimed < 0` or the claimant is not in `tubule_meta`; one not in `tubule_meta` sets it where `claimed < 0`. `counts_by_compartment`, a dict of `all`, `intra` and `extra`, each of shape `(staged tubules, genes, N_BINS)`, from one `np.bincount(gene * N_BINS + bin)` per tubule and compartment, where `intra` is the transcripts whose `cell_id` is not -1 and `tube_gene_bin = counts_by_compartment["all"]`, with `bin = np.clip(np.searchsorted(bin_edges, position, side="left") - 1, 0, N_BINS - 1)`. `stage_by_tubule` is `stage_of.loc[tubule_order]` |
| G4 | `GENE_SIGMA = (np.floor(n_tubules / 100 * 1.5), 1.5)`; `total_by_bin = tube_gene_bin.sum(axis=1)`; `transcript_density = total_by_bin / (tubule_area[:, None] / N_BINS) * 1e3`; `equal_total = total_by_bin * (TUBE_LIBRARY_SIZE / total_by_bin.sum(axis=1))[:, None]`, `TUBE_LIBRARY_SIZE = 1e6`; `gaussian_filter(arr, sigma=GENE_SIGMA, mode=["wrap", "reflect"])`. `plot_total_transcripts` draws the density, the smoothed density and the smoothed equal total, `coolwarm`, each on its own scale |
| G5 | `equal_genes = tube_gene_bin * (TUBE_LIBRARY_SIZE / per_tubule)[:, None, None]`, his `normalize_tube_counts_dict`; `MARKER_GENES` maps a layer to a panel gene; each gene's `(tubules, bins)` array is smoothed as in G4, his `smooth_tube_gene_bin`. `plot_gene_heatmaps` draws unsmoothed above smoothed, `coolwarm`, each panel on its own scale |
| G6 | `library_factor = TUBE_LIBRARY_SIZE / per_tubule`, from all transcripts; for each compartment `smoothed_by_compartment[key]` is each gene's `smooth(counts[:, gene, :] * library_factor[:, None])`. `drawn` is the genes with at least `MIN_GENE_TRANSCRIPTS` (500) transcripts in the staged tubules. `across_radius[key] = smoothed.mean(axis=0)[drawn]`, genes by bins; `radius_order` sorts the genes by their centre of mass across the bins in the `all` panel. `plot_all_genes` draws the three compartments with that order, each row z-scored, `coolwarm`, `vmin=-3`, `vmax=3`, and the markers labelled |
| G7 | `across_cycle[key] = smoothed.mean(axis=2).T[drawn]`, genes by tubules; `ring_position` is the angle of each gene's circular centre of mass over `ring_angle = 2 * pi * arange(n_tubules) / n_tubules`, `arctan2(sum(w * sin), sum(w * cos)) % (2 * pi)` on the `all` panel, and `cycle_order` sorts by it. Drawn with `plot_all_genes` and the stage bar above |

**Phase H**

| Step | Code |
|---|---|
| H1 | `staged_cells` is the cells whose `tubule_id` is in `stage_of`. `composition(column, names, cells)` is `pd.crosstab(cells["tubule_id"], cells[column])` divided by its row sums, once for `cell_type` and once for `germ_type` among the germ cells; the mean per stage is `shares.groupby(stage_of.reindex(shares.index)).mean()`; `scipy.stats.kruskal` per type across the four groups. Every staged tubule is used; the `good_cells_bin` filter is not applied |
| H2 | The staged tubule per stage whose `gc_shares` row has the smallest summed absolute difference from the stage mean; `read_crop` reads its padded bounding box from the finest image level, DAPI and Cellbound3 only, and `show_tubule` fills its cells by type |

Later analyses are meant to read each cell's fine type and drop cells where `excluded` is set. Only the first half exists: `germ_type` is the fine type of a germ cell, and nothing is excluded, because the cluster steps that would set it are deferred. In Phase E, per-cell
results (`tubule_id`, `in_tubule`) go in the table. The tubule polygons are loaded in the
notebook's E1 as `tubules_um`, and Phase G reads them from there.

## Not built yet

Rows and blocks planned in `PLAN.md` under *Where to go deeper*. Numbered as they will be.

| Step | Code |
|---|---|
| C7 | *Not run; the notebook does not cluster.* For each resolution `r`: `sc.tl.leiden(adata, resolution=r, key_added=f"leiden_{r}", flavor="igraph", n_iterations=2, random_state=0)`, then our own `relabel_by_size()` renumbers the clusters 1..k by size |
| D1 | *Deferred.* `cluster_key = "leiden_1.0"` |
| D2 | `sc.tl.rank_genes_groups(adata, groupby=cluster_key, layer="lognorm", method="wilcoxon")`, read back with `sc.get.rank_genes_groups_df` |
| D7 | *Deferred.* `pd.crosstab(adata.obs[cluster_key], adata.obs["cell_type"], normalize="index")`, kept in `adata.uns`; `germ_share = composition[GERM].sum(axis=1)`, where `GERM = ["SPG", "Scytes", "STids", "Elongating"]` |
| D8 | *Deferred.* `cluster_type` is `"Germ"` where `germ_share` is at least `PURITY` (0.90), even when one germ type alone is; otherwise the top type where its share is at least `PURITY`, and `"Ambiguous"` where neither is, written to `adata.obs["cluster_type"]`; how it was made to `adata.obs["cluster_call_source"]`, one of `auto`, `manual` or `unresolved`. `cell_type` is left as the per-cell call, with all four germ types |
| D9 | *Deferred.* A manual call is written to `cluster_type` with source `manual`. A cluster left unresolved keeps `Ambiguous` with source `unresolved`, and its cells get `adata.obs["excluded"] = True` and `adata.obs["excluded_at"] = "major"` |
| E4 | *Not built.* `tubules` is `tubule_meta` plus `n_cells` and `germ_share` per tubule from `adata.obs.groupby("tubule_id")`, with columns `good_cells`, `good_circle`, `good_cells_bin`, `area`, `eccen`, `maj_axis_length`. Areas and lengths are in mosaic pixels |
| E5 | *Not built.* `plot_tubules_on_tissue` and `plot_tubule_qc`; the share inside by cell type is `adata.obs.groupby("cell_type", observed=True)["in_tubule"].mean()` |
| H3 | *Not built.* `counts` is the H1 crosstab without normalizing; `area_um2 = staged["area"] * UM_PER_PX ** 2` (`UM_PER_PX = 0.108`); `density = counts.div(area_um2, axis=0) * 1e4`, cells per 10,000 µm², with `germ` and `all cells` columns added. Needs the `staged` table of E4 and F1 (manual) |

**Radial layering (rest of it)**

| Step | Code |
|---|---|
| R3 | Share of cells with `inside` false, per tubule and overall (11 of 131,377 on b2r0); share of rays whose `hit` is not a single LineString, and the tubules they come from (not yet counted); distance from `polygon.centroid` to the median of the tubule's cell centres (median 10 µm on b2r0, 0.07 of the equivalent radius). A star-shape check is not needed for the position itself; if wanted, rebuild the polygon from the outline's points sorted by angle, and compare its area with the outline's (98.6% of b2r0's tubules agree within 2%) |
| R4 | `adata.obs.groupby("cell_type", observed=True)["radial_pos"]` and the same by GC cluster, drawn as violins or boxes, ordered by the expected order |
| R5 | *Remaining:* genes by their full tubule-by-bin pattern, clustered; genes with two peaks; and the comparison with his `*-all-counts-dict.pkl` (not staged). The markers, the total and every gene averaged over one axis are G3 to G7 |

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
