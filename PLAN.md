# Analysis plan

Dated 2026-11-02. Each of the collaborator's 16 mouse testis sections is analysed on its own: QC,
clustering and cell-type annotation, then tubules identified, layered and staged. That lets
cell-type proportions, gene and signalling dynamics, and the interstitium be compared by stage,
and finally between young and old animals. Phases A to H run per section;
only I pools sections.

The analysis is worked out end to end on one section, b2r0, before it touches the others:
phases A to H are settled there first, then run on the remaining 15 sections, and I follows
once every section is through H. Anything learned on the later sections that changes a step
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

### Notebook sections for Phases B and C

The notebook's first sections, in order, with the two checkpoints marked:

1. Load the section's table
2. Filter cells
3. Normalize
4. PCA, neighbours and UMAP
5. Leiden clustering at a few resolutions: 0.5, 1.0, 1.5 and 2.0

   **Checkpoint 1, clustered:** the filtered table with a cluster column per resolution.

6. Choose the resolution to annotate, 1.0 to start
7. Marker genes per cluster, at that resolution
8. Per-cell correlation against each reference
9. Cluster composition of the per-cell calls
10. Cluster calls

    **Checkpoint 2, annotated:** the whole table, with everything sections 1 to 10 added: the
    clusters at every resolution, the markers, every `corr_<type>` score, the per-cell calls and
    the cluster calls. The later phases read whichever columns they need from it, so which
    label each one uses is decided when that phase is built.

The split lets clustering and annotation be iterated on separately. Changing the resolution,
the markers, the correlation or the calls reruns only sections 6 to 10 from checkpoint 1;
changing the clustering reruns from section 2 and invalidates checkpoint 2 with it, since its
key includes checkpoint 1's. The chosen resolution is part of checkpoint 2's key, not
checkpoint 1's, so trying another of the four never reclusters. The cluster calls follow from
the composition and the purity threshold, so the threshold is part of checkpoint 2's key too.

Four resolutions rather than one because extra Leiden runs on the same neighbour graph are
cheap, and because 1.0 may be too coarse for the Green reference: the earlier b2r0 report had
19 clusters at 1.0, several Ambiguous, against Green's 32 types, 16 of them germ cells. Germ
cells are a continuous trajectory that Leiden cuts at arbitrary points, so what a resolution is
judged by is whether its clusters are fine enough to match the reference's types.

**One object passes through the whole notebook.** Section 1 loads the section's table as a
single AnnData, and every later section adds to that same object rather than making a new one:
raw counts kept in `layers["counts"]` before normalization overwrites `X`, embeddings in
`obsm`, clusters, scores and calls in `obs`, markers and parameters in `uns`. A checkpoint is
that object written to disk, and loading one puts the notebook back where it was. The
SpatialData store's images and boundaries are not loaded into it; tissue figures read the
shapes they need from the store when they draw.

Memory clean-up is not planned. Adding to one object in place leaves no stale copies to free,
and the object is small: about 147,000 cells by 491 genes is under 300 MB dense, plus the
neighbour graph. Two habits keep it that way: no `adata.copy()` per section, and anything big
that is only needed briefly, such as a full reference dataset, is deleted once used.

A third checkpoint could go after normalization (section 3), so that changing PCA, neighbours
or Leiden does not redo filtering and normalization. Not needed yet; add it if clustering
b2r0 turns out slow.

**The notebook starts on CPU (scanpy).** b2r0 has about 147,000 cells after filtering, over
491 genes, which CPU clusters in minutes, and a CPU job runs on any OSCER node or locally with
little or no wait.

GPU (rapids-singlecell, as step 2 used) is the option if clustering becomes the slow part. It
clusters faster, but OSCER's GPU queue is the wait: step 2's jobs on `sooner_gpu_test` were
scheduled hours apart, the last a day out. The `wsl` machine's own GPU has no queue.
rapids-singlecell mirrors scanpy's API, so moving over is a switch in the notebook rather than
a rewrite. The two Leiden implementations do not give identical clusters, so the backend is
part of checkpoint 1's key.

## Phase A — Segmentation and data QC

**Scope:** Per section.

**Deferred.** The segmentations already run are accepted as they are, and the analysis moves on.
The report's existing QC section stays as a first look. Nothing is filtered on QC yet beyond
Phase B's low count cut.

### Sweeps for later

Not now. QC checks that could become thresholds, and how to judge where to set them. Each
applies to the segmentation as it stands; none means re-segmenting.

| QC check | Default | Could try | How to judge it |
|---|---|---|---|
| Cell volume | No cut | Drop the smallest and largest cells, by fixed volumes or quantiles | Whether small cells are fragments (few genes, no clear type) and large ones merges (markers of mutually exclusive types) |
| Genes per cell | No cut | A minimum, alongside Phase B's count cut | Whether the cells it drops cluster together or are scattered through every type |
| Blank codewords per cell | No cut | A maximum share of a cell's counts | Blanks against real genes per section; whether high-blank cells sit in one area of tissue |
| Transcripts assigned to a cell | Not checked | A minimum share per section or per FOV | Share per FOV drawn on tissue; a low area points at a segmentation or imaging fault |
| Cells with no polygon | Kept | Dropped | How many there are (`has_boundary`); whether their calls differ from cells with one |
| Boundary plane | Every plane kept | Drop cells whose boundary is on the top or bottom plane (`obs["z_plane"]`) | Whether those cells have fewer counts or different calls |
| FOV seams and tissue edges | Kept | Drop cells within a set distance of a seam or edge | The QC metrics drawn on tissue; whether seam and edge cells are lower in counts or odd in type |
| Damaged areas | Kept | Mask regions by hand or by low transcript density | The same tissue maps |
| Outlier sections | All kept | Flag a section whose metrics sit apart from the rest | The metrics compared across sections, pooled |

Every row is judged by the same two questions: are the cells it removes worse than the rest,
and does removing them change the clusters or the Phase C calls.

## Phase B — Filtering and clustering

**Scope:** Per section, with one default recipe for every section.

**The priority is to reach cell-type annotations (Phase C).** Phases B and C run on defaults,
untuned, so the analysis can move on to the tubule phases. Where a sweep could improve them, and
how to judge it, is noted at the end of each section and left for later.

The notebook starts from the recipe step 2 (`bin/cluster_spatialdata_gpu.py`) uses, as a
default for every section, clustered at four resolutions and annotated at 1.0.

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
| Leiden | Resolutions 0.5, 1.0, 1.5 and 2.0, all written; 1.0 is the one annotated; random seed 0 throughout |

Each resolution leaves two columns. `_v0` is Leiden's own labelling; `_v1` renumbers clusters
1..k by size, so cluster 1 is the largest. **v1 is what every call and figure names a cluster
by.** The report adds a render-only `v2`, the cluster's position down the dendrogram, so
neighbouring numbers are transcriptionally similar; figures read `v1 (v2)`.

Normalizing to each section's own median makes log-normalized values a per-section scale.
Anything compared across sections starts from counts.

Step 4 (`bin/create_centroids.py`) sums counts and CP10K per cluster at every resolution.
Sums pool, so any union of clusters is the sum of its rows, and `log1p(X / n_cells)` is
comparable to the `ln(mean + 1)` reference centroids.

Marker genes per cluster are not part of this phase: they are computed at the resolution being
annotated, as the first section of Phase C's half of the notebook.

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
| Resolution | 0.5, 1.0, 1.5, 2.0; 1.0 annotated | 0.1 to 2.0 in steps of 0.1 | Agreement between neighbouring resolutions (adjusted Rand index); whether a split gives clusters with distinct markers; how pure each cluster's Phase C calls are |
| Seed | 0 | Several | Agreement across seeds at the chosen resolution |

Two checks apply to most rows: how stable the clusters are, and how cleanly Phase C can name them.
A cluster that only appears at one setting, or that no reference type takes, is the thing to
look for. Whether thresholds and resolution end up fixed for the cohort or set per section is
decided by the same sweeps.

## Phase C — Cell-type annotation

**Scope:** Per section, with the same references and the same purity threshold for the whole
cohort.

Every cell is correlated against a reference, and each cluster is then named from the mix of
its cells' calls alone. The notebook follows step 3's method (`bin/annotate_celltypes.py`) for
the per-cell calls, on defaults.

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
5. **Cluster call.** A cluster is named when one type takes at least 90% of its cells. Below
   that it is called Ambiguous, and gets its own slide laying out the evidence. The earlier
   report used 70% (`PURITY_THRESHOLD`); 90% asks for a much clearer majority.
6. **Cluster markers, alongside.** Each cluster's top marker genes, from notebook section 7,
   are shown next to its composition to see whether they are informative. They are not used to
   make or check a call, and no marker list is defined for the reference cell types.

| Reference | Species | Cell types | Values | Fit for this plan |
|---|---|---|---|---|
| Green et al. 2018 | Mouse | 11 major types; 32 in the collaborator's copy: 7 somatic, 9 Sertoli subtypes, 4 SPG, 12 germ-cell clusters | ln(mean TP10K + 1); 19 of the collaborator's 32 columns confirmed against GEO | Used first. Resolves the germ-cell stages that layering and staging need |
| Shami et al. (`shami_human_testis_centroids.csv.gz`) | Human | 11 major types | ln(mean + 1) | Run alongside. Too coarse for staging; useful as a cross-check |

**Both references are used, Green first.** Shami is run alongside as a cross-check.

Green et al. 2018, *A Comprehensive Roadmap of Murine Spermatogenesis Defined by Single-Cell
RNA-Seq*, Dev Cell 46(5):651-667 ([PubMed 30146481](https://pubmed.ncbi.nlm.nih.gov/30146481/)):
34,633 Drop-seq mouse testis cells. Its data is GEO
[GSE112393](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112393), checked on
2026-11-03:

| File | Size | Holds |
|---|---|---|
| `GSE112393_MergedAdultMouseST25_DGE.txt.gz` | 81 MB | Raw UMI counts, 37,241 genes by 34,633 cells |
| `GSE112393_MergedAdultMouseST25_PerCellAttributes.txt.gz` | 1.5 MB | Per cell: dataset, `CellType` as one of the 11 major types, and QC metrics |
| `GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx` | 11 MB | The authors' own centroids for germ-cell clusters GC1-12, over 24,475 genes, as average, total and fraction expressing |

The 11 major types are InnateLymph, Macrophage, Endothelial, Myoid, Leydig, Sertoli, Unknown,
SPG, Scytes, STids and Elongating. The paper normalizes as `ln(TP10K + 1)`: counts divided by
the cell's total UMIs, times 10,000, then natural log. A centroid is `ln(mean TP10K + 1)`, the
same convention as Shami's.

**What the collaborator's 32 columns turn out to be:**

| Columns | Source | Status |
|---|---|---|
| 7 somatic (InnateLymph, Macrophage, Endothelial, Myoid, Leydig, Sertoli, Tcf21+) | The major types, recomputed from the counts and labels above | Exact match on every panel gene; Tcf21+ is the GEO type called Unknown |
| `Germ_1` to `Germ_12` | GC1-GC12 in the spreadsheet, rounded to 3 decimals | Exact match |
| `SPG_1` to `SPG_4`, 9 Sertoli subtypes | Not on GEO: per-cell labels below the 11 major types were not deposited | Unconfirmed; source unknown |

So the scale question is settled for 19 of 32 columns, and they need no recomputing: the
collaborator's values are the published ones. The 13 subtype columns cannot be rebuilt from
public data.

**The reference is built once, outside the per-section notebook.** A separate notebook downloads
these files and writes centroid tables to `assets/reference/`, the same shape as Shami's, over
the 491 panel genes. Every per-section notebook then loads those small tables rather than the
full dataset. Two levels come out of it:

- **Green, 11 major types,** recomputed from the counts and labels. Coarse, and every column is
  reproducible from public data.
- **Green, fine,** the 7 somatic types and GC1-12 as published. Whether the 13 unconfirmed SPG
  and Sertoli subtype columns join it is still to decide.

Seven panel genes are missing from the counts under the panel's names, and nine from the
spreadsheet, being newer symbols for renamed genes (`Prxl2a`, `Plpp3`, `Mesd`, `Meioc`, ...).
The reference notebook maps them to their older names rather than dropping them.

A possible extra, not planned yet: reference marker genes computed from the reference data
itself, rather than a list written by hand. How depends on what the reference comes as. With
per-cell counts and labels, markers come from a differential expression test between types;
with centroids alone, from how specific each gene is to one type relative to the rest. Either
way they are limited to the 491 panel genes, and could be set beside each cluster's own
markers to see whether the two agree.

Calls must mean the same thing in every section for Phase I. The same references and the same
90% threshold in every section are what make them comparable.

### Sweeps for later

Not now. Where the defaults could be swept, and what would show a choice is better:

| Parameter | Default | Could try | How to judge it |
|---|---|---|---|
| Reference | Green and Shami, both run | Green's 32 types grouped to coarser levels | Agreement on the coarse types both have, once Green's types are grouped to Shami's 11 |
| Correlation | Spearman | Pearson on log counts, as the collaborator's script does | Share of cells whose call changes, and which types they move between |
| Genes | All genes the panel and reference share | Only the genes that vary most across the reference's types, or the computed reference markers | Agreement with the all-gene calls |
| Low-confidence cells | Every cell gets its top type | Leave a cell uncalled when its top two types are close | Share of cells left uncalled; whether those are low-count or on boundaries |
| Purity threshold | 90% | 70% to 95% | How many clusters become Ambiguous, and whether their own markers set them apart |

Normalization needs no sweep here: scaling or log-transforming a cell's counts keeps their
ranks, so Spearman gives the same answer on any of them.

Whatever the setting, a call can be sanity-checked by where the type sits in tissue, such as
germ cells inside tubules and Leydig cells outside them.

## Phase D — Identify tubules

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

## Phase E — Radial layering

**Scope:** Per section.

Each cell in a tubule gets a radial position: 0 at the basement membrane, 1 at the lumen.

- Expected order, outside in: spermatogonia and Sertoli nuclei, spermatocytes, round
  spermatids, elongating spermatids.
- Checks Phases C and D together; a wrong order points at one of them.
- Gene expression along the radial axis, within each cell type.

## Phase F — Stage the tubules

**Scope:** Per section. Could also be pooled: ordering tubules around the cycle may be steadier
with every section's tubules together than with one section's alone.

Each tubule cross-section sits at one stage of the epithelial cycle (I-XII in mouse), set by
which germ cells occur together.

- **Rule-based:** germ-cell associations and stage markers (e.g. Stra8 in stages VII-VIII).
- **Unsupervised:** a composition or pseudobulk profile per tubule, ordered around the cycle;
  tubule profiles often form a ring.
- Realistic target: stage groups (I-III, IV-VI, VII-VIII, IX-XII) rather than all twelve.
- Ask the collaborator whether these sections have been staged before; that is a ready validation set.

## Phase G — Within the tubule, by stage

**Scope:** Per section.

- Cell-type proportions by stage
- Gene expression across the cycle per cell type, as pseudobulk per tubule ordered by stage
- Sertoli-germ cell signalling by stage, from ligand-receptor pairs

Only pairs on the 491-gene panel can be tested; listing which ones are is a cheap early step.

## Phase H — Extratubular cells and the tubule

**Scope:** Per section.

- Interstitial types: Leydig, macrophages (peritubular vs interstitial), endothelial, peritubular
  myoid
- Neighbourhood enrichment, and each cell's distance to the nearest tubule
- Whether the interstitium beside a tubule changes with that tubule's stage
- Ligand-receptor signalling across the basement membrane

## Phase I — Young vs old

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
- [ ] Green reference: whether the unconfirmed SPG and Sertoli subtype columns are used, and if so where they came from
- [ ] Which image channels exist for tubule segmentation?
- [ ] Where tubule polygons are stored, given the store is write-once
- [ ] Has the collaborator staged any of these sections?
- [ ] Which ligand-receptor pairs the panel covers
