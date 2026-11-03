# Analysis plan

## The data

Each sample is one section of mouse testis, measured with spatial transcriptomics: a panel of
491 genes read out in place across the tissue. Cells were segmented in 3D with cellpose by the
lab's own pipeline, and each transcript was assigned to the cell it falls in. There are 16
sections, from young and old animals. The analysis below is written for one section; only the last
phase brings sections together.

## The biology

Sperm are made inside the seminiferous tubules, and a cross-section of mouse testis cuts through
many of them. Inside each tubule, germ cells mature from the outside in. Spermatogonia sit
against the basement membrane; spermatocytes above them pass through meiosis; round spermatids
and then elongating spermatids lie towards the lumen, into which sperm are released. Sertoli
cells span that whole depth and support every stage. Outside the tubules, peritubular myoid
cells wrap each one, and the interstitium between them holds Leydig cells, which make
testosterone, along with macrophages and blood vessels.

The epithelium does not mature uniformly. It runs through a repeating cycle of twelve stages,
I to XII in mouse, and each stage is a fixed combination of germ cells at different steps of
development. A single tubule cross-section sits at one stage, so one section of testis is a
sample of the whole cycle, spread across its tubules.

## The analysis

Because expression is measured in every cell in place, a section can be read at all of
these levels at once. For a section, the analysis:

1. **Identifies every cell** by its expression: which somatic or germ cell type it is (Phases B
   to D).
2. **Finds the tubules,** places each cell inside or outside one, and orders the cells within a
   tubule from the basement membrane to the lumen, which should recover the germ cells'
   developmental order (Phases E and F).
3. **Stages each tubule** from the germ cells it contains (Phase G).
4. **Asks how the tubule changes through the cycle:** which cell types are present at each stage,
   how their gene expression shifts, and how Sertoli and germ cells signal to each other
   (Phase H).
5. **Asks how the interstitium relates to the tubule beside it:** where Leydig cells and
   macrophages sit, and whether their neighbourhood or their signalling across the basement
   membrane follows that tubule's stage (Phase I).

Repeated over sections from young and old animals, the same measurements show how ageing
changes the cycle, the tubule and the interstitium (Phase J).

Each phase that is decided lists its steps in order, and what each one decides. `NOTEBOOK.md`
is how the notebook that runs them is built, with the code and the column names for each step.

## Phase A — Build the stores

**Scope:** Per section.

Each section's raw output and cellpose segmentation are turned into a SpatialData store (images,
transcripts and cell boundaries) and a table beside it (counts and per-cell metadata). This is
the pipeline's steps 1a and 1b, run before the notebook; `CLAUDE.md` has the detail. Every later
phase reads the table, and the store only to draw tissue.

A section is ready for Phase B once its table names the shapes in its store. Two sections need
extra care: c1r1's segmentation raster is cropped, and a1r0 holds two tissue sections in one
region, so it is split first.

## Phase B — Segmentation and data QC

**Scope:** Per section.

**Minimal.** The segmentations already run are accepted as they are. This phase only looks: a few
plots per section, so a bad region or an odd distribution is seen before clustering. Nothing is
filtered here; the only cut is Phase C's low count cut.

### Plots

1. **Transcripts per cell and genes per cell,** as histograms.
2. **Cell volume,** as a histogram.
3. **Blank codewords against real genes,** as the share of each cell's counts that are blanks.
4. **Transcript count drawn on tissue,** to show FOV seams, edges and damaged areas at a glance.

Anything that looks wrong is noted for the sweeps below, not acted on in this phase.

### Sweeps for later

Not now. QC checks that could become thresholds, and how to judge where to set them. Each
applies to the segmentation as it stands; none means re-segmenting.

| QC check | Default | Could try | How to judge it |
|---|---|---|---|
| Cell volume | No cut | Drop the smallest and largest cells, by fixed volumes or quantiles | Whether small cells are fragments (few genes, no clear type) and large ones merges (markers of mutually exclusive types) |
| Genes per cell | No cut | A minimum, alongside Phase C's count cut | Whether the cells it drops cluster together or are scattered through every type |
| Blank codewords per cell | No cut | A maximum share of a cell's counts | Blanks against real genes per section; whether high-blank cells sit in one area of tissue |
| Transcripts assigned to a cell | Not checked | A minimum share per section or per FOV | Share per FOV drawn on tissue; a low area points at a segmentation or imaging fault |
| Cells with no polygon | Kept | Dropped | How many there are; whether their calls differ from cells with one |
| Boundary plane | Every plane kept | Drop cells whose boundary is on the top or bottom plane | Whether those cells have fewer counts or different calls |
| FOV seams and tissue edges | Kept | Drop cells within a set distance of a seam or edge | The QC metrics drawn on tissue; whether seam and edge cells are lower in counts or odd in type |
| Damaged areas | Kept | Mask regions by hand or by low transcript density | The same tissue maps |
| Outlier sections | All kept | Flag a section whose metrics sit apart from the rest | The metrics compared across sections, pooled |

Every row is judged by the same two questions: are the cells it removes worse than the rest,
and does removing them change the clusters or the Phase D calls.

## Phase C — Filtering and clustering

**Scope:** Per section, with one default recipe for every section.

**The priority is to reach cell-type annotations (Phase D),** so Phases C and D run on defaults,
untuned. Where a sweep could improve them, and how to judge it, is noted at the end of each
phase and left for later.

### Steps

1. **Load the section's table,** keeping the raw counts alongside the normalized values.
2. **Filter cells.** Drop cells under 20 transcripts. The high cut, for doublets and merged
   cells, is off; if turned on it is a quantile of transcript count, taken before the low cut.
   All 491 genes are kept: no highly-variable selection and no gene filter, since the panel is
   already curated.
3. **Normalize.** Scale every cell to the section's median total, then log. Keep the
   log-normalized values for markers.
4. **Scale.** Unit variance, not centred, clipped at 10.
5. **PCA,** 30 components.
   - Plot: variance explained per component, to show 30 PCs are enough.
6. **Neighbour graph,** 15 neighbours, cosine distance.
7. **UMAP.**
   - Plot: UMAP coloured by transcript count, to show whether structure tracks depth rather
     than type.
8. **Leiden at 0.5, 1.0, 1.5 and 2.0.** Each resolution's clusters are renumbered 1..k by size,
   so cluster 1 is the largest. Every call and figure names a cluster by that number.
   - Plot: UMAP coloured by cluster, one panel per resolution.
   - Plot: cells per cluster at each resolution.
   - Plot: the clusters drawn on tissue at the resolution annotated, to show whether they are
     spatially coherent.

The same recipe, steps 3 to 8, is used again when a single cell type is subclustered in
Phase D step 9, restarted from its raw counts rather than reusing the section-wide result.

1.0 may turn out too coarse once annotation narrows down to finer types: germ cells are a
continuous trajectory that Leiden cuts at arbitrary points, so a resolution is judged by whether
its clusters are fine enough to match the reference's types.

### Sweeps for later

Not now. Where the defaults could be swept, and what would show a choice is better:

| Parameter | Default | Could try | How to judge it |
|---|---|---|---|
| Low count cut | 20 transcripts | 10 to 50 | Whether any cluster is mostly low-count cells; whether the dropped cells map to one area of tissue or are scattered |
| High count cut | Off | Top 0.5% or 1% | Whether the cells it removes co-express markers of mutually exclusive types, as doublets and merged cells would |
| Normalization | Section median, then log1p | CP10K; dividing by cell volume | Whether clusters still track cell size or transcript count rather than type |
| PCs | 30 | 20 to 50 | Variance explained per PC; how much the clusters change |
| Neighbours | k = 15 | 10 to 50 | How much the clusters change |
| Resolution | 0.5, 1.0, 1.5, 2.0; 1.0 annotated | 0.1 to 2.0 in steps of 0.1 | Agreement between neighbouring resolutions (adjusted Rand index); whether a split gives clusters with distinct markers; how pure each cluster's Phase D calls are |
| Seed | 0 | Several | Agreement across seeds at the chosen resolution |

Two checks apply to most rows: how stable the clusters are, and how cleanly Phase D can name them.
A cluster that only appears at one setting, or that no reference type takes, is the thing to
look for. Whether thresholds and resolution end up fixed for the cohort or set per section is
decided by the same sweeps.

## Phase D — Cell-type annotation

**Scope:** Per section, with the same reference and the same purity threshold for the whole
cohort, which is what makes calls comparable across sections for Phase J.

Every cell is correlated against a reference, and each cluster is then named from the mix of
its cells' calls alone.

### Steps

0. **Build or load the reference, once for the cohort.** The first section builds it and the
   rest load it. Green's 11 major types are rebuilt from the published counts and labels, as the
   mean expression per type over the panel genes. Green's GC1-12 are taken as published, and
   Shami's table is loaded as it stands. Panel genes renamed since a reference was made are
   mapped to the names it uses rather than dropped. The reference used is recorded with the
   annotation.
1. **Choose the resolution,** 1.0 to start.
2. **Marker genes per cluster,** by a Wilcoxon test on the log-normalized values.
   - Plot: a dot plot of each cluster's top markers.
3. **Per-cell correlation.** Spearman correlation of each cell's raw counts against every
   reference type, over the genes the panel and reference share. Genes are matched
   case-insensitively, so mouse `Acta2` finds human `ACTA2`.
   - Plot: the number of genes shared between the panel and the reference.
4. **Standardize within each cell.** A cell's correlation to every type rises with the genes it
   captured; standardizing makes cells comparable with each other without changing any cell's
   top type.
5. **Per-cell call.** The type each cell correlates with most.
6. **Cluster composition.** Each cluster's share of cells per called type.
   - Plot: a heatmap of that composition, with the purity threshold marked.
7. **Cluster call.** A cluster is named for its top type when that type takes at least 90% of
   its cells, and Ambiguous otherwise. Each call records how it was made: automatically, by
   hand, or left unresolved.
   - Plot: UMAP coloured by cluster call, with Ambiguous clusters picked out.
8. **Resolve Ambiguous clusters by hand.** Each Ambiguous cluster is investigated and given a
   type if the evidence supports one. The evidence is its composition from step 6, its top
   markers from step 2, and where its cells sit in tissue. A short note of the evidence is kept
   per cluster. One rule for what counts as enough evidence applies to every section, so manual
   calls mean the same thing across the cohort. A cluster that stays unresolved keeps
   Ambiguous, and its cells are excluded at the major level. Every later phase drops excluded
   cells. The excluded share per section is reported, and compared across sections in Phase J.
   - Plot: the calls drawn on tissue, to check where each type sits, such as germ cells inside
     tubules and Leydig cells outside them.
   - Plot: a bar chart of cells per call and of how many were excluded.
9. **Subcluster the types that have a finer reference.** Once the major types are settled, take
   the cells of one major type and run the Phase C recipe on them again, then annotate the
   subclusters against that type's finer reference by steps 3 to 8: the same correlation, the
   same composition rule, and the same by-hand resolution, recorded the same way. Which types
   have a finer reference:

   | Major type | Fine reference | Status |
   |---|---|---|
   | Germ cells (SPG, Scytes, STids, Elongating) | Green GC1-12 | **Now** |
   | Sertoli, SPG subtypes | The collaborator's file, if cleared for use | Open question |
   | Leydig, Macrophage, Endothelial, Myoid, InnateLymph, Unknown (Tcf21+) | None yet | Later; other references to find |

   The result is each cell's fine cell type, **which is what every later phase reads.** A type
   with no fine reference keeps its major label, and cells in subclusters left unresolved by
   hand are excluded.
   - Plots: steps 2 to 8 are drawn again at the fine level.

### References

**One reference at a time.** The reference is a setting, Green or Shami, and a run annotates
against that one only. Comparing the annotations the two give is left for later; each
reference's annotation is kept, so that comparison needs no rerun.

| Reference | Species | Cell types | Values | Use |
|---|---|---|---|---|
| Green et al. 2018, 11 major types | Mouse | InnateLymph, Macrophage, Endothelial, Myoid, Leydig, Sertoli, Unknown (Tcf21+), SPG, Scytes, STids, Elongating | ln(mean TP10K + 1), rebuilt from the published counts and labels | **The first runs.** |
| Green et al. 2018, germ-cell clusters | Mouse | GC1-GC12 | ln(mean TP10K + 1), as published | Narrowing down germ cells (D9), for staging |
| Shami et al. | Human | 11 major types | ln(mean + 1) | The alternative |

Green et al. 2018, *A Comprehensive Roadmap of Murine Spermatogenesis Defined by Single-Cell
RNA-Seq*, Dev Cell 46(5):651-667 ([PubMed 30146481](https://pubmed.ncbi.nlm.nih.gov/30146481/));
data at GEO [GSE112393](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112393). Each
level is a reference of its own, never mixed into one table.

### Sweeps for later

Not now. Where the defaults could be swept, and what would show a choice is better:

| Parameter | Default | Could try | How to judge it |
|---|---|---|---|
| Reference | Green's 11 major types | Shami | Agreement between the two annotations, cell by cell and cluster by cluster |
| Reference level | Major level, then narrowing down (step 9) | Green's fine level for every cell at once | Whether the fine calls agree with the narrowed-down ones |
| Correlation | Spearman | Pearson on log counts, as the collaborator's script does | Share of cells whose call changes, and which types they move between |
| Genes | All genes the panel and reference share | Only the genes that vary most across the reference's types, or reference markers computed from the reference data | Agreement with the all-gene calls |
| Low-confidence cells | Every cell gets its top type | Leave a cell uncalled when its top two types are close | Share of cells left uncalled; whether those are low-count or on boundaries |
| Purity threshold, major and fine level | 90% | 70% to 95%, each level on its own | How many clusters become Ambiguous; the share of all cells, and of germ cells, left excluded once hand resolution is done; whether an Ambiguous cluster's own markers set it apart |

Whatever the setting, a call can be sanity-checked by where the type sits in tissue, such as
germ cells inside tubules and Leydig cells outside them.

## Phases E to J

These are an outline. Phases A to D are decided; E to J are planned in detail once the major-type
annotation is running. Every one of them reads each cell's fine cell type and drops excluded
cells.

## Phase E — Identify tubules

**Scope:** Per section.

Every cell gets a tubule and a compartment (tubule, boundary or interstitium), and every tubule
a polygon.

- **From cells:** peritubular myoid cells (Acta2, Myh11) trace each boundary; germ and Sertoli
  cells fill the inside. Connected tubular cells on a spatial graph, closed with a concave hull.
- **From the image:** segment tubules and lumens from the mosaic (DAPI, polyT or a boundary
  stain).
- Likely a combination: cells for the boundary, the image to check it.
- QC per tubule: area, circularity, cell count. Flag tubules cut obliquely, since they distort
  radial position and staging.

Per-cell results are stored with the cells; tubule polygons are not per cell, so they are
stored separately.

## Phase F — Radial layering

**Scope:** Per section.

Each cell in a tubule gets a radial position: 0 at the basement membrane, 1 at the lumen.

- Expected order, outside in: spermatogonia and Sertoli nuclei, spermatocytes, round
  spermatids, elongating spermatids.
- Checks Phases D and E together; a wrong order points at one of them.
- Gene expression along the radial axis, within each cell type.

## Phase G — Stage the tubules

**Scope:** Per section. Could also be pooled: ordering tubules around the cycle may be steadier
with every section's tubules together than with one section's alone.

Each tubule cross-section sits at one stage of the epithelial cycle (I-XII in mouse), set by
which germ cells occur together.

- **Rule-based:** germ-cell associations and stage markers (e.g. Stra8 in stages VII-VIII).
- **Unsupervised:** a composition or pseudobulk profile per tubule, ordered around the cycle;
  tubule profiles often form a ring.
- GC1-12 are germ-cell differentiation states, not cycle stages. Staging uses them as input, as
  the germ-cell composition of each tubule, and does not read a stage from them directly.
- Realistic target: stage groups (I-III, IV-VI, VII-VIII, IX-XII) rather than all twelve.

## Phase H — Within the tubule, by stage

**Scope:** Per section.

- Cell-type proportions by stage
- Gene expression across the cycle per cell type, as pseudobulk per tubule ordered by stage
- Sertoli-germ cell signalling by stage, from ligand-receptor pairs on the 491-gene panel

## Phase I — Extratubular cells and the tubule

**Scope:** Per section.

- Interstitial types: Leydig, macrophages (peritubular vs interstitial), endothelial, peritubular
  myoid
- Neighbourhood enrichment, and each cell's distance to the nearest tubule
- Whether the interstitium beside a tubule changes with that tubule's stage
- Ligand-receptor signalling across the basement membrane

## Phase J — Young vs old

**Scope:** Over the cohort, sections grouped by animal and animals by age.

- Distribution of tubules across stages
- Composition within each stage
- Tubule size
- Interstitium: Leydig and macrophage numbers and position
- Expression by stage and cell type

The unit of replication is the **animal**, not the section or the tubule: pseudobulk per animal,
or mixed models with animal as a random effect. Sections from one animal are not independent.
Normalizing to each section's own median makes log-normalized values a per-section scale, so
anything compared across sections starts from counts.

## Decisions and open questions

- [ ] Which sections are young and which old, and how many animals in each? (Do the a1, b2, c1,
      d2, 2a, 2b prefixes mark animals?)
- [ ] Green SPG and Sertoli subtypes: the collaborator's file has them, but they are not in the public data. Use them, and if so, where did they come from?
- [ ] Which image channels exist for tubule segmentation?
- [ ] Has the collaborator staged any of these sections? That would be a ready validation set for Phase G.
- [ ] Which ligand-receptor pairs the panel covers
