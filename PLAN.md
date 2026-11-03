# Analysis plan

Dated 2026-11-02. Each of the collaborator's 16 mouse testis sections is analysed on its own: QC,
clustering and cell-type annotation, then tubules identified, layered and staged. That lets
cell-type proportions, gene and signalling dynamics, and the interstitium be compared by stage,
and finally between young and old animals. Phases A to F run per section;
only G pools sections.

The analysis is worked out end to end on one section, b2r0, before it touches the others:
phases A to F are settled there first, then run on the remaining 15 sections, and G follows
once every section is through F. Anything learned on the later sections that changes a step
goes back to b2r0 too, so all 16 are analysed the same way.

The per-section analysis is one Quarto notebook rather than the `bin/` scripts, built up a
section at a time in phase order. It starts from the store and table steps 1a and 1b already
wrote; the steps after those are a reference for the recipe, not what runs.

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

`NOTES.md` tracks where the runs stand; this file is what the analysis is for and in what order.

## Phase A1 — Segmentation and data QC

**Scope:** Per section.

**Deferred.** The segmentations already run are accepted as they are, and the analysis moves on.
The report's existing QC section stays as a first look. Ideas for when this is picked up:

- Cells per section; transcripts, genes and volume per cell
- Share of transcripts assigned to a cell; blank codewords against real genes
- Which z-plane each boundary came from (`obs["z_plane"]`), and cells with no polygon
  (`has_boundary`)
- The same metrics drawn on tissue, to catch FOV seams, edge effects and damaged areas
- The metrics compared across sections, pooled, to show an outlier section

## Phase A2 — Filtering and clustering

**Scope:** Per section, with one default recipe for every section.

**The priority is to reach cell-type annotations (A3).** A2 and A3 run on defaults, untuned,
so the analysis can move on to the tubule phases. Where a sweep could improve them, and how to
judge it, is noted at the end of each section and left for later.

The notebook starts from the recipe step 2 (`bin/cluster_spatialdata_gpu.py`) uses, as a
default for every section, read at resolution 1.0.

| Stage | Default |
|---|---|
| Filter, low | Drop cells under 20 transcripts (`--min_counts`) |
| Filter, high | Off (`--max_counts_quantile`); when on, the top quantile of transcript counts is dropped for doublets and merged cells, taken before the low cut |
| Genes | All panel genes kept: no highly-variable selection and no gene filter, since the 491-gene panel is already curated |
| Normalize | Total counts scaled to the section's median cell total, then log1p |
| Scale | Unit variance, not centred, clipped at 10 |
| PCA | 30 components |
| Neighbors | Cosine distance, rapids default k = 15 |
| UMAP | On that graph |
| Leiden | Resolutions 0.1 to 2.0 in steps of 0.1, all written; 1.0 is the one read; random seed 0 throughout |

Each resolution leaves two columns. `_v0` is Leiden's own labelling; `_v1` renumbers clusters
1..k by size, so cluster 1 is the largest. **v1 is what every call and figure names a cluster
by.** The report adds a render-only `v2`, the cluster's position down the dendrogram, so
neighbouring numbers are transcriptionally similar; figures read `v1 (v2)`.

Normalizing to each section's own median makes log-normalized values a per-section scale.
Anything compared across sections starts from counts.

Step 4 (`bin/create_centroids.py`) then sums counts and CP10K per cluster at every resolution.
Sums pool, so any union of clusters is the sum of its rows, and `log1p(X / n_cells)` is
comparable to the `ln(mean + 1)` reference centroids.

Still to do now:

- Add top marker genes per cluster; the report does not show them yet.

### Sweeps for later

Not now. Where the defaults could be swept, and what would show a choice is better:

| Parameter | Default | Could try | How to judge it |
|---|---|---|---|
| Low count cut | 20 transcripts | 10 to 50 | Whether any cluster is mostly low-count cells; whether the dropped cells map to one area of tissue or are scattered |
| High count cut | Off | Top 0.5% or 1% | Whether the cells it removes co-express markers of mutually exclusive types, as doublets and merged cells would |
| Other cell filters | None | Cell volume, genes per cell | The same checks as the count cuts |
| Normalization | Section median, then log1p | CP10K; dividing by cell volume | Whether clusters still track cell size or transcript count rather than type |
| PCs | 30 | 20 to 50 | Variance explained per PC; how much the clusters change |
| Neighbours | k = 15 | 10 to 50 | How much the clusters change |
| Resolution | 1.0 | 0.1 to 2.0, already written | Agreement between neighbouring resolutions (adjusted Rand index); whether a split gives clusters with distinct markers; how pure each cluster's A3 calls are |
| Seed | 0 | Several | Agreement across seeds at the chosen resolution |

Two checks apply to most rows: how stable the clusters are, and how cleanly A3 can name them.
A cluster that only appears at one setting, or that no reference type takes, is the thing to
look for. Whether thresholds and resolution end up fixed for the cohort or set per section is
decided by the same sweeps.

## Phase A3 — Cell-type annotation

**Scope:** Per section, against one reference and one rulebook fixed for the whole cohort.

Every cell is correlated against a reference, and each cluster is then named by a person from
the mix of its cells' calls. The notebook follows step 3's method (`bin/annotate_celltypes.py`),
on defaults, which makes no cluster-level call itself.

1. **Per-cell correlation.** Spearman correlation of each
   cell's raw counts against every reference cell type, over the genes the panel and reference
   share. Genes match case-insensitively, so mouse `Acta2` finds human `ACTA2`; the matched genes
   are written to `<sample>.gene_overlap.tsv`.
2. **Standardize within each cell.** A cell's correlation to every type rises with the genes it
   captured, so correlations are z-scored per cell. This cannot change a cell's top type; it
   makes cells comparable with each other.
3. **Per-cell call.** `cell_type_per_cell` is the type with the highest correlation; every
   type's value is kept as `corr_<type>`.
4. **Cluster composition (report).** Each cluster's share of cells per called type.
5. **Cluster call (by hand).** A cluster is named when one type takes at least 70% of its cells
   (`PURITY_THRESHOLD`). Below that it is called Ambiguous, and gets its own slide laying out
   the evidence.

| Reference | Species | Cell types | Values | Fit for this plan |
|---|---|---|---|---|
| Green et al. 2018, from the collaborator (`green_mouse_testis_centroids.csv.gz`) | Mouse | 32: 7 somatic, 9 Sertoli subtypes, 16 germ-cell types | Unconfirmed | Resolves the germ-cell stages that layering and staging need |
| Shami et al. (`shami_human_testis_centroids.csv.gz`, current default) | Human | 11 major types | ln(mean + 1) | Too coarse for staging; useful as a cross-check |

**Recommendation:** annotate with the Green mouse reference, once the collaborator confirms how
its values were computed. Their own script uses Pearson on log2 counts per 1,000; ours is
Spearman, which uses only ranks, so the scale matters less but still needs confirming.

Calls must mean the same thing in every section for Phase G. One fixed reference does most of
that; the hand-made cluster calls follow one written rulebook (marker genes per type, the 70%
threshold).

### Sweeps for later

Not now. Where the defaults could be swept, and what would show a choice is better:

| Parameter | Default | Could try | How to judge it |
|---|---|---|---|
| Reference | Green (mouse) | Shami (human) | Agreement on the coarse types both have, once Green's 32 types are grouped to Shami's 11 |
| Correlation | Spearman | Pearson on log counts, as the collaborator's script does | Share of cells whose call changes, and which types they move between |
| Genes | All genes the panel and reference share | Only each type's reference markers | Agreement with the all-gene calls; marker expression in the cells called |
| Low-confidence cells | Every cell gets its top type | Leave a cell uncalled when its top two types are close | Share of cells left uncalled; whether those are low-count or on boundaries |
| Purity threshold | 70% | 50% to 90% | How many clusters become Ambiguous, and whether their markers say they should be |

Normalization needs no sweep here: scaling or log-transforming a cell's counts keeps their
ranks, so Spearman gives the same answer on any of them.

Whatever the setting, a call is checked the same two ways: canonical markers should be
expressed in the cells called that type, and the type should sit where it belongs in tissue,
such as germ cells inside tubules and Leydig cells outside them.

## Phase B — Identify tubules

**Scope:** Per section.

Every cell gets a `tubule_id` and a `compartment` (tubule, boundary or interstitium), and every
tubule a polygon.

- **From cells:** peritubular myoid cells (Acta2, Myh11) trace each boundary; germ and Sertoli
  cells fill the inside. Connected tubular cells on a spatial graph, closed with a concave hull.
- **From the image:** segment tubules and lumens from the mosaic (DAPI, polyT or a boundary
  stain).
- Likely a combination: cells for the boundary, the image to check it.
- QC per tubule: area, circularity, cell count. Flag tubules cut obliquely, since they distort
  radial position and staging.

Where the polygons live is open: the store is written once by design, so a tubule shapes element
goes beside the table (GeoParquet, say) or needs an exception. Per-cell columns fit in the table.

Prototype in a notebook on b2r0; promote to a pipeline step once the method settles.

## Phase C — Radial layering

**Scope:** Per section.

Each cell in a tubule gets a radial position: 0 at the basement membrane, 1 at the lumen.

- Expected order, outside in: spermatogonia and Sertoli nuclei, spermatocytes, round
  spermatids, elongating spermatids.
- Checks Phases A3 and B together; a wrong order points at one of them.
- Gene expression along the radial axis, within each cell type.

## Phase D — Stage the tubules

**Scope:** Per section. Could also be pooled: ordering tubules around the cycle may be steadier
with every section's tubules together than with one section's alone.

Each tubule cross-section sits at one stage of the epithelial cycle (I-XII in mouse), set by
which germ cells occur together.

- **Rule-based:** germ-cell associations and stage markers (e.g. Stra8 in stages VII-VIII).
- **Unsupervised:** a composition or pseudobulk profile per tubule, ordered around the cycle;
  tubule profiles often form a ring.
- Realistic target: stage groups (I-III, IV-VI, VII-VIII, IX-XII) rather than all twelve.
- Ask the collaborator whether these sections have been staged before; that is a ready validation set.

## Phase E — Within the tubule, by stage

**Scope:** Per section.

- Cell-type proportions by stage
- Gene expression across the cycle per cell type, as pseudobulk per tubule ordered by stage
- Sertoli-germ cell signalling by stage, from ligand-receptor pairs

Only pairs on the 491-gene panel can be tested; listing which ones are is a cheap early step.

## Phase F — Extratubular cells and the tubule

**Scope:** Per section.

- Interstitial types: Leydig, macrophages (peritubular vs interstitial), endothelial, peritubular
  myoid
- Neighbourhood enrichment, and each cell's distance to the nearest tubule
- Whether the interstitium beside a tubule changes with that tubule's stage
- Ligand-receptor signalling across the basement membrane

## Phase G — Young vs old

**Scope:** Over the cohort, sections grouped by animal and animals by age.

- Distribution of tubules across stages
- Composition within each stage
- Tubule size
- Interstitium: Leydig and macrophage numbers and position
- Expression by stage and cell type

The unit of replication is the **animal**, not the section or the tubule: pseudobulk per animal,
or mixed models with animal as a random effect. Sections from one animal are not independent.

## Decisions and open questions

- [ ] Which sections are young and which old, and how many animals in each? (Do the a1, b2, c1,
      d2, 2a, 2b prefixes mark animals?)
- [ ] Reference: confirm with the collaborator how the Green centroids were computed, then annotate with it
- [ ] Which image channels exist for tubule segmentation?
- [ ] Where tubule polygons are stored, given the store is write-once
- [ ] Has the collaborator staged any of these sections?
- [ ] Which ligand-receptor pairs the panel covers
