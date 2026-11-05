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
2. **Finds the tubules** and places each cell inside or outside one (Phase E).
3. **Stages each tubule** from the germ cells it contains (Phase F).
4. **Places each cell between the middle of its tubule and the wall,** and shows where each cell
   type sits there across the cycle (Phase G).
5. **Asks how the tubule changes through the cycle:** which cell types are present at each stage
   (Phase H).

Planned and not yet built, under *Where to go deeper*: QC plots, checks on the radial coordinate, the gene heatmaps, how expression and signalling shift through the
cycle, how the interstitium relates to the tubule beside it, and the same measurements across
young and old animals.

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

### Steps

1. **Load the section's table,** keeping the raw counts, and count each cell's transcripts over
   the real genes. Blank codewords are not counted.
2. **The tissue image.** The DAPI stain from the section's store, read at a coarse level and
   placed on the cells' micrometre axes, as the ground truth every later tissue figure is read
   against.
Anything that looks wrong is noted for the sweeps below, not acted on in this phase.

### Sweeps for later

Not now. QC checks that could become thresholds, and how to judge where to set them. Each
applies to the segmentation as it stands; none means re-segmenting.

| QC check | Default | Could try | How to judge it |
|---|---|---|---|
| Cell volume | No cut | Drop the smallest and largest cells, by fixed volumes or quantiles | Whether small cells are fragments (few genes, no clear type) and large ones merges (markers of mutually exclusive types) |
| Genes per cell | No cut | A minimum, alongside Phase C's count cut | Whether the cells it drops cluster together or are scattered through every type |
| Blank codewords per cell | No cut | A maximum share of a cell's counts | The blank counts per cell from Phase B step 5; whether high-blank cells sit in one area of tissue, from the step 6 map |
| Transcripts assigned to a cell | Not checked | A minimum share per section or per FOV | Share per FOV drawn on tissue; a low area points at a segmentation or imaging fault |
| Cells with no polygon | Kept | Dropped | How many there are; whether their calls differ from cells with one |
| Boundary plane | Every plane kept | Drop cells whose boundary is on the top or bottom plane | Whether those cells have fewer counts or different calls |
| FOV seams and tissue edges | Kept | Drop cells within a set distance of a seam or edge | The QC metrics drawn on tissue; whether seam and edge cells are lower in counts or odd in type |
| Damaged areas | Kept | Mask regions by hand or by low transcript density | The same tissue maps |
| Outlier sections | All kept | Flag a section whose metrics sit apart from the rest | The metrics compared across sections, pooled |

Every row is judged by the same two questions: are the cells it removes worse than the rest,
and does removing them change the clusters or the Phase D calls.

## Phase C — Filtering and dimension reduction

**Scope:** Per section, with one default recipe for every section.

**The priority is to reach cell-type annotations (Phase D),** so Phases C and D run on defaults,
untuned. Where a sweep could improve them, and how to judge it, is noted at the end of each
phase and left for later.

### Steps

1. **Filter cells.** Drop cells under 20 transcripts. The high cut, for doublets and merged
   cells, is off; if turned on it is a quantile of transcript count, taken before the low cut.
   All 491 genes are kept: no highly-variable selection and no gene filter, since the panel is
   already curated.
2. **Normalize.** Scale every cell to the section's median total, then log. Keep the
   log-normalized values for markers.
3. **Scale.** Unit variance, not centred, clipped at 10.
4. **PCA,** 30 components.
   - Plot: variance explained per component, to show 30 PCs are enough.
5. **Neighbour graph,** 15 neighbours, cosine distance.
6. **UMAP.**
   - Plot: UMAP coloured by transcript count, to show whether structure tracks depth rather
     than type.
The notebook draws neither plot yet; the UMAP is drawn in Phase D, coloured by the calls.

Nothing is clustered. Cells are annotated one at a time in Phase D, and the UMAP is for looking
at the data. If clusters are needed later, for the deferred cluster steps in Phase D or for
reclustering a cell type, they are Leiden at 0.5, 1.0, 1.5 and 2.0, each
resolution's clusters renumbered 1..k by size so that cluster 1 is the largest, and a cluster is
always named by that number. The same recipe, steps 2 to 6 then Leiden, would be restarted from
raw counts on a single cell type's cells.

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

Every cell is correlated against a reference, and its top type is its call.

**The first pass annotates cells, not clusters.** Each cell is called germ or somatic from its
own call. Nothing is clustered, and nothing is named at the cluster level. On b2r0 the per-cell
calls separate the two cleanly on the UMAP, and drawn on tissue they put germ cells inside the
tubules and somatic cells between them. The per-cell call cannot cleanly separate neighbouring
germ stages, so a second level takes the germ cells alone and correlates them against Green's
GC1-12 (step 6). That level is correlation again on a finer reference, not a reclustering.

### Steps

Steps are numbered as in the notebook. D1, D2 and D7 to D9 are cluster steps that are not run, so
the notebook goes from D0 to D3 and stops at D6. They are under *Where to go deeper*.

0. **Load the reference, made once for the cohort.** The table is made outside the notebook and
   every section loads it; the notebook never builds it, and a missing table is an error.
   Green's 11 major types are rebuilt from the published counts and labels, as the
   mean expression per type over the panel genes. Green's GC1-12 are taken as published, and
   Shami's table is loaded as it stands. Panel genes renamed since a reference was made are
   mapped to the names it uses rather than dropped. The reference used is recorded with the
   annotation.
3. **Per-cell correlation.** Spearman correlation of each cell's raw counts against every
   reference type, over the genes the panel and reference share. Genes are matched
   case-insensitively, so mouse `Acta2` finds human `ACTA2`.
   - Print: the number of genes shared between the panel and the reference.
4. **Standardize within each cell.** A cell's correlation to every type rises with the genes it
   captured; standardizing makes cells comparable with each other without changing any cell's
   top type.
5. **Per-cell call.** The type each cell correlates with most, and from it germ or somatic: a
   cell is germ when its call is SPG, Scytes, STids or Elongating, and somatic otherwise. This
   is the first pass's annotation.
   - Plot: UMAP and tissue coloured by the per-cell type, then by germ or somatic.
6. **Germ-cell subtypes.** The cells called germ, correlated against Green's GC1-12 by steps 3
   to 5 again, with nothing reclustered. GC1-12 are differentiation states in order, so on
   tissue they should run from the wall inward. A cell that is not germ has no germ type. The
   result is each germ cell's fine type, `germ_type`, and every other cell keeps its major type.
   - Plot: UMAP and tissue coloured by germ-cell type.
### References

**One reference at a time.** The reference is a setting, Green or Shami, and a run annotates
against that one only. Comparing the annotations the two give is left for later; each
reference's annotation is kept, so that comparison needs no rerun.

| Reference | Species | Cell types | Values | Use |
|---|---|---|---|---|
| Green et al. 2018, 11 major types | Mouse | InnateLymph, Macrophage, Endothelial, Myoid, Leydig, Sertoli, Unknown (Tcf21+), SPG, Scytes, STids, Elongating | ln(mean TP10K + 1), rebuilt from the published counts and labels | **The first runs.** |
| Green et al. 2018, germ-cell clusters | Mouse | GC1-GC12 | ln(mean TP10K + 1), as published | Germ-cell subtypes (D6), for staging |
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
| Reference level | Major level, then GC1-12 on germ cells (step 6) | Green's fine level for every cell at once | Whether the fine calls agree with the two-level ones |
| Correlation | Spearman | Pearson on log counts, as the collaborator's script does | Share of cells whose call changes, and which types they move between |
| Genes | All genes the panel and reference share | Only the genes that vary most across the reference's types, or reference markers computed from the reference data | Agreement with the all-gene calls |
| Low-confidence cells | Every cell gets its top type | Leave a cell uncalled when its top two types are close | Share of cells left uncalled; whether those are low-count or on boundaries |
| Germ grouping | Four germ types pooled to name a cluster | Each type on its own, the earlier rule | Share of cells in clusters that pass: on b2r0 at 1.0 it is 26% alone and 84% pooled |
| Purity threshold, major and fine level | 90% | 70% to 95%, each level on its own | How many clusters become Ambiguous; the share of all cells, and of germ cells, left excluded once hand resolution is done; whether an Ambiguous cluster's own markers set it apart |

Whatever the setting, a call can be sanity-checked by where the type sits in tissue, such as
germ cells inside tubules and Leydig cells outside them.

## Phase E — Identify tubules

**Scope:** Per section.

Every cell gets a tubule, or none, and every tubule carries a few quality flags.

**This phase uses Alex Vargo's existing results rather than a new segmentation.** His tubule
boundaries are traced from the image: the dark interior of each tubule is found, touching
tubules are cut apart by a morphological opening, and the wall of each is followed along the
bright Cellbound3 ring (with DAPI) by a minimum-cost path. He then assigned every cell to the
tubule whose boundary contains its centre, and wrote a per-tubule metadata table. We take his
cell-to-tubule map and metadata as inputs for every section. Rebuilding his method is later work.

A cell-only method was tried first on b2r0 and abandoned: a distance graph and a Delaunay graph
on the germ cells both merge touching tubules into pieces of tens of thousands of cells, and a
density watershed seeds the wrong places because tubules are rings with large lumens.

### Steps

Numbered as in the notebook.

1. **Load the outlines and the flags.** One `<id>.npz` per tubule, moved into micrometres, and
   each tubule's flags from his metadata: a quality flag for the cells, and whether the outline is
   a good circle. Each tubule gets one of four flags.
2. **Draw the tubules on the section.** The outlines over DAPI, coloured by flag, and over the
   per-cell type and the germ-cell type.
3. **Load the map.** The map has one entry per cell of the section's table, in the table's order,
   so its length must equal the table's; a cell outside every tubule has none. Each cell gets a
   tubule id from it.
The tubule table and the check of the assignment are not built; they are under *Where to go
deeper*.

Not done here:

- **Compartments.** Only inside or outside a tubule. A boundary compartment needs the polygons.
- **Oblique tubules.** His eccentricity and circle flags are carried through; a rule for
  excluding obliquely cut tubules is left until radial layering shows whether they distort the radial
  order.
- **The polygons.** Radial layering needs them for the radial position. They are one `<id>.npz` per
  tubule, holding a `global_path` array of `(row, column)` points in mosaic pixels. The notebook
  already loads them and moves them into micrometres (E1 in the notebook).

## Phase F — Stage the tubules

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

**For now, Alex's stages are used.** He staged tubules by eye from per-tubule images, into
exactly those four groups (a few sections also have finer values such as 3-4 or 7), and the
calls are in the tubule metadata of five sections: a1r0 upper and lower, b2r2, c1r1 and d2r0.
The other eleven, including b2r0, have none.

**His computed stage covers every section.** `all_tubule_order_250113.csv` places each good-circle
tubule of all sixteen sections on a ring, from a principal component analysis of the tubules'
radial gene profiles pooled across sections. It gives `ang`, the angle on the ring, `order`, the
rank around it pooled across sections, and `prop_stage`, a group from 1 to 4 cut from the angle.
His cutoffs were fit on his by-eye calls, which they match 78 to 92% exactly and almost always
within one group. We load these and do not recompute them; b2r0 has 77 staged tubules (17, 28,
11 and 21 in groups 1 to 4). `ang` is circular, so any statistic on it must respect the wrap. For
plots along x, `order` is re-ranked within the section. The phases below run on the staged
sections, and a staging of our own is later work to be scored against his.

### Steps

Numbered as in the notebook.

1. **Load each tubule's stage and its position on the ring.** His group, 1 to 4, and the
   continuous `phase` from his angle, with the four cuts between groups.
2. **Draw the tubules by stage,** by group and by position on the ring, on DAPI.
3. **Draw the stage on the UMAP and the tissue,** each cell taking its tubule's group and then its
   tubule's position on the ring.

**Showing the continuous stage.** The tubule outlines on the section are drawn once by group and
once by `ang`, so the two can be read against each other.
- `ang` is rebased to start at the beginning of group 1 and run 0 to 2π, and coloured with a
  cyclic colormap, because the cycle wraps: group 4 sits next to group 1.
- The colourbar carries the four group cuts, taken from the smallest and largest rebased `ang` in
  each group.
- The same colour is given to the cells of each tubule on the tissue and the UMAP, as the group
  colour already is.
- Read it for gradients: neighbouring tubules at neighbouring stages. A group whose tubules sit
  at both ends of the bar, or a jump inside one group, shows where his cuts fall in the
  continuous variable.

A first look on b2r2 (140 tubules) shows how far germ composition alone gets. From each tubule's
shares of the 12 GC clusters, a regularized logistic regression recovers his group with 76%
cross-validated accuracy, against 49% for always guessing the biggest group. Late stages are
clean (9-12: 95% of tubules right) and 7-8 is the weak one (32%; 63% with balanced class
weights, at the cost of 4-6 falling from 83% to 53%). Nearly every error lands on the adjacent
group, so composition recovers the order of stages and blurs the boundaries between his groups.
These are exploratory numbers from one section, not a method.

## Phase G — Radial position

**Scope:** Per section. Every cell in a tubule gets a position, staged or not; the heatmaps draw the
staged tubules.

Each cell gets a radial position, from the middle of its tubule out to its wall, by Alex Vargo's
definition (`radialTubule.create_tubule_transcript_polar_info`). **The direction is his: 0 at the
centroid, 1 at the wall.** It uses his outlines, in micrometres from Phase E.

### Steps

Numbered as in the notebook.

1. **Radial position of every cell.** `rad` is the distance from the outline's centroid to the cell.
   `dist_to_edge` is found by casting a ray from the cell, outward along its angle from the centroid,
   for the length of the outline's bounding-box diagonal: the length of the part of that ray inside the
   outline, and where the ray leaves and re-enters, the distance to the farthest crossing. A cell not
   strictly inside the outline, so one on its boundary too, gets 0. `radial_pos` is
   `rad / (rad + dist_to_edge)`. The ray is cast from the cell, so nothing assumes the outline is
   star-shaped. On b2r0, 131,377 cells are in a tubule, and 11 are not strictly inside their own outline.
2. **Cell-type density by radial position.** The cells of each type counted into 50 bins of equal area,
   with edges at the square root of k over 50, per staged tubule. The density is the count over the bin's
   area, the outline's area over 50, in cells per 10,000 µm². It is not rescaled to the same total for each
   tubule, so a more crowded tubule stays more crowded.
   - Plot: one heatmap per type, for SPG, Scytes, STids, Elongating, Sertoli and Leydig. Tubules run
     along x in their order on the stage ring, with the stage groups marked above, and the radial bin runs
     up from the middle at the bottom. The top row is the density as counted, the bottom row smoothed.
   - The smoothing is a Gaussian, σ = 2 tubules and 1.5 bins, wrapping across tubules because the ring
     closes and reflecting across bins. σ is set by eye for 77 tubules; his b2r0 notebooks use 1 to 2.
     Each panel runs from 0 to its own 99th percentile.

3. **Transcripts in each tubule, by radial position.** The section's transcripts, blank codewords
   dropped. Each goes in the tubule whose outline contains it, boundary included, and gets the same
   radial position as a cell. A transcript in more than one outline is kept by the tubule his rule gives
   it (`xscript-to-tubule-assignment.py`): tubules in order of their id; a tubule gets every transcript in
   its outline that nobody has claimed; a tubule he flagged good for cells also takes transcripts claimed
   earlier; a tubule in his table but not good takes only those claimed by a tubule not in his table; a
   tubule not in his table takes only unclaimed ones. Counts are kept as tubules by genes by bins, for the
   staged tubules, in three compartments: all transcripts, those inside a cell, and those outside any cell
   (30.1% on b2r0), by the cell the store gives each transcript. On b2r0, 90.3% of the 34.8 million transcripts of the panel's genes are in a tubule, and
   12.1 million of them in the 77 staged ones.
4. **Total transcripts by radial position.**
   - Plot: all genes together, tubules along x in ring order and radial bins up from the middle. Two
     panels are a density, transcripts per 1,000 µm² of the bin, unsmoothed and smoothed; the third is his
     normalization, every tubule rescaled to the same total, then smoothed.
   - The smoothing is his: a Gaussian whose width across tubules is `floor(n_tubules / 100 * 1.5)`, so 1
     for 77 tubules, and 1.5 across bins, wrapping across tubules and reflecting across bins.
5. **Marker genes by radial position.** His gene heatmaps, rebuilt from our own transcripts: each tubule
   rescaled to the same total, each gene smoothed as in step 4, one marker for each layer, each panel on
   its own colour scale. The markers are the panel genes most specific to their type in Green's
   reference: Stra8 (SPG), Ccna1 (Scytes), Acrv1 (STids), Prm1 (Elongating) and Amhr2 (Sertoli).

6. **Every gene across the radius.** To look for spatial patterns not already known. Each gene is one row:
   its counts across the radial bins, averaged over the staged tubules, normalized and smoothed as in step 5.
   The inside and the outside compartments are scaled by the tubule's all-transcript factor, as in his
   notebook, so they add up to the total. Each row is on its own z-score scale, and the genes are sorted by
   where their counts sit across the radius, in the order of the all-transcripts panel. Genes with fewer than
   500 transcripts in the staged tubules are left out: 467 of 491 on b2r0.
   - Plot: three panels, all transcripts, inside a cell and outside any cell, with the same gene order.
7. **Every gene around the cycle.** The same for the ring: each gene's counts across the tubules in ring
   order, averaged over the radius, sorted by the angle at which they sit on the ring.
   - The ring has no natural start, so a gene high at both ends of the plot, like the spermatocyte genes, is
     one peak that wraps round. A gene with two separate peaks sits between them in this order, so it can look
     like a gene with one broad peak. The marginal also hides how the radius and the cycle combine.

What it shows on b2r0: Sertoli and SPG at the wall in every tubule, with Leydig at the very edge, as
cells inside generous outlines; STids at mid-radius through the early part of the cycle; Elongating towards
the lumen in the later part; and Scytes at mid-to-outer radius, dense at both ends of the cycle. Median
`radial_pos` for good circles runs Elongating 0.55, STids 0.64, Scytes 0.79, SPG 0.92, Sertoli 0.93.

Total transcripts are densest at mid-radius, a median of about 3,400 per 1,000 µm² in bins 20 to 29,
against about 2,000 toward the lumen and 2,300 at the wall, and the lumen end is emptiest in tubules about
25 to 45 round the ring. The markers fall where their layer sits: Stra8 near the wall in one stage group,
Ccna1 at mid-radius at both ends of the cycle, Acrv1 at mid-radius through the middle of it, Prm1 toward the lumen
in a few tubules, and Amhr2 at the wall in every tubule. Our per-tubule totals match his
`b2r0-all-counts-dict.pkl` to within 0.1% for all 77 staged tubules, and the radial profile of each
correlates with his at 0.999 or better. The Elongating cell-type density peaks in a different part of the
cycle from Prm1, which has not been looked into.

Across the radius the genes form a clean diagonal, from lumen genes (Prm1) through mid-radius (Acrv1, Ccna1) to
wall genes (Amhr2, Stra8), a continuous gradient and not sharp layers. It is tightest for transcripts inside
cells. Outside cells the order is the same but noisier, with an edge band in the outermost bins that many wall genes
light up, and the lumen genes also show up outside cells in the innermost bins. Around the cycle the genes
form a diagonal wave, each peaking at a ring position, clearest inside cells; the median gene's range around the
cycle is about 1.3 times its mean overall, 1.4 inside cells and 1.65 outside.

Known limits:

- **No lumen.** Position is scaled to the centroid, not to the lumen, so the cells nearest the lumen do not
  reach 0 and a layer sits at a different value in a tubule with a large lumen.
- **Tubules that are not star-shaped.** A ray that crosses the wall more than once is measured to its
  farthest crossing, and the equal-area bins are exact only for star-shaped outlines. On b2r0, 405 of
  131,366 rays (0.31%) cross more than once, in 39 tubules and nearly all unstaged.
- **Raw bins are sparse.** A bin holds about 1,100 µm² and a few cells.

## Phase H — Composition by stage

**Scope:** Per section, on sections that have stage calls. Developed on b2r2, which is fully
staged; b2r0 is not.

### Steps

Numbered as in the notebook.

1. **Composition by stage.** For each staged tubule, each cell type's share of its cells, and
   each germ-cell type's share of its germ cells, averaged by stage group, with a Kruskal-Wallis
   test per type. His staging is independent of our calls, so this is not circular even for germ
   cells. The notebook uses every staged tubule; the flagged-good filter on the cells is not
   applied yet.
   - Plot: stacked bars of the mean composition per stage group, by cell type and by germ-cell
     type.
2. **One tubule per stage.** The staged tubule whose germ-cell shares are closest to its stage's
   mean, drawn at full resolution with its cells filled by type over DAPI and Cellbound3. They
   illustrate the means and do not test them.
What b2r2 shows (125 tubules: 14 in 1-3, 61 in 4-6, 18 in 7-8, 32 in 9-12): the germ types move
with stage as expected, with STids falling from 35-39% of cells to 6%, Elongating peaking at 45% in
7-8, and Scytes high at both ends of the cycle. Sertoli rise from 7.1% of cells in 4-6 to 11.1% in
9-12 and their density rises with them, from 8.4 to 11.8 per 10,000 µm², so that change is not just
the denominator. SPG rise weakly. Total density falls in the late groups, 117 to 119 in groups
1-6 down to 104 in 9-12. The other interstitial types inside the outlines do not vary with stage
and carry no information about it: they are cells that fall inside generous boundaries, not
residents of the tubule.

The density figures in that paragraph come from the exploration and are not in the notebook. Caveats: one section; 14 and 18 tubules in the smaller groups; the 4-6 group is wide, and its
cell densities spread widely, as if it held several stages; the outline area includes the lumen, so
density mixes cell packing with lumen size; p values are not corrected for the number of types
tested.

## Where to go deeper

Everything above is built and runs in the notebook. What follows is planned and not built, in the
order the work would go. Each item says what it adds and what it depends on. Everything here reads
each cell's fine type, and would drop excluded cells once the cluster steps below set them.

### QC plots (Phase B)

3. **Transcripts per cell and genes per cell,** as histograms.
4. **Cell volume,** as a histogram.
5. **Blank codewords against real genes,** in two panels: each codeword's total counts across the
   section, blanks beside the real genes, to show how far the genes sit above the noise floor; and
   the number of blank counts each cell collects.
6. **Transcript count drawn on tissue,** to show FOV seams, edges and damaged areas at a glance. The
   mean transcripts per cell along x and along y is drawn at the map's margins, and the imaging
   fields of view are marked on all three, so a dip can be matched to a field-of-view boundary or
   told apart from a smooth fall towards the tissue's edge.

### Tubule table and assignment check (Phase E)

4. **Tubule table.** One row per tubule: its cells, the share of them
   that are germ, and his metadata (area, eccentricity and major axis length).
5. **Check the assignment.** *Not in the notebook yet.*
   - Plot: the tubules on tissue, each in its own colour, with germ cells that fall outside any
     tubule picked out.
   - Plot: cells per tubule, and the share of each cell type that falls inside a tubule. Germ
     and Sertoli cells should mostly be inside, Leydig and the other interstitial types mostly
     outside.

What the check shows on b2r0: 143 tubules of a median 700 cells; 92% of germ cells and 87% of
Sertoli cells inside a tubule. The boundaries are generous, so about half of Leydig and
Endothelial cells also land inside one. The germ cells left outside, 7.6%, are mostly the
low-count strip along one edge of the section and a few tubules that have no outline.

### Radial layering (rest of it)

The coordinate and the cell-type density heatmap are Phase G. What is left, numbered as in the earlier
plan:

3. **Check the outlines and the coordinate.**
   - Report the share of cells that are not strictly inside their own outline, which get position 1.
   - Report the share of rays that cross the wall more than once, and the tubules they come from.
   - Report the distance between each outline's centroid and the median position of its cells.
4. **Check the order of the cell types.** Outside in, the expected order is Sertoli and SPG at the
   wall, then spermatocytes, round spermatids and elongating spermatids, with the interstitial
   types at the very edge. The GC clusters should run roughly in index order from the wall inward.
   - Plot: radial position by cell type, and by GC cluster.
   - This checks Phases D and E together. A wrong order points at one of them.
5. **The rest of the gene heatmaps.** Markers and every gene, averaged over one axis, are drawn in Phase G.
   What is left:
   - Genes by their full tubule-by-bin pattern, not one axis at a time, clustered to find domains, as his
     `lap_labs` do.
   - Genes with two peaks, which the sorted marginals hide.
   - The check against his `*-all-counts-dict.pkl` is made in the kernel and is not in the notebook, because the
     pickle is not staged.

`radial_pos` is the coordinate used. It does not depend on how large a tubule is. His code also
keeps `rad` and `dist_to_edge`, the distances to the centroid and along the ray to the wall, in
micrometres, and they are kept for the same reason: a cell type's `dist_to_edge` grows with the
tubule's size, while `radial_pos` does not.

The limits are listed under Phase G. A lumen boundary would sharpen the Sertoli/SPG and spermatid
layers; expression along the radial axis within a cell type is a later step that reads `radial_pos`.

### More by stage (Phase H)

3. **Cell density by stage.** Cells of each type per 10,000 µm² of the
   tubule's outline, to separate a real change in a cell type from a change in the share of the
   others.
   - Plot: boxplots per stage group for the four germ types, Sertoli and all cells.

Also:

- Gene expression across the cycle per cell type, as pseudobulk per tubule ordered by stage
- Sertoli-germ cell signalling by stage, from ligand-receptor pairs on the 491-gene panel

### The interstitium and the tubule

**Scope:** Per section.

- Interstitial types: Leydig, macrophages (peritubular vs interstitial), endothelial, peritubular
  myoid
- Neighbourhood enrichment, and each cell's distance to the nearest tubule
- Whether the interstitium beside a tubule changes with that tubule's stage
- Ligand-receptor signalling across the basement membrane

### Young against old

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

### Cluster steps (Phase D)

Not run, since the first pass annotates cells. They name clusters from the mix of their cells'
calls and need a clustering that Phase C does not run. Numbered as in the notebook.

1. **Choose the resolution,** 1.0 to start. Needed only by the cluster steps 2 and 7 to 9.
2. **Marker genes per cluster,** by a Wilcoxon test on the log-normalized values. *Deferred.*
   - Plot: a dot plot of each cluster's top markers.
7. **Cluster composition.** Each cluster's share of cells per called type, and its
   share of germ cells, the four germ types (SPG, Scytes, STids, Elongating) taken together. The
   per-cell calls keep all four, and are pooled only to name a cluster: a cluster between two
   germ stages splits its cells between them, whatever its resolution.
   - Plot: a heatmap of that composition, with the purity threshold marked twice: a box where
     one type reaches it, and a box of another colour across the germ columns where the germ
     cells together do.
8. **Cluster call.** A cluster is named Germ when the germ types together take at
   least 90% of its cells, including a cluster that is almost all one germ type, and otherwise
   for its top type when that type takes at least 90%. Anything else is Ambiguous. Each call
   records how it was made: automatically, by hand, or left unresolved.
   - Plot: UMAP coloured by cluster call, with Ambiguous clusters picked out.
9. **Resolve Ambiguous clusters by hand.** Each Ambiguous cluster is investigated and
   given a type if the evidence supports one. The evidence is its composition from step 7, its
   top markers from step 2, and where its cells sit in tissue. A short note of the evidence is
   kept per cluster. One rule for what counts as enough evidence applies to every section, so
   manual calls mean the same thing across the cohort. A cluster that stays unresolved keeps
   Ambiguous, and its cells are excluded at the major level. Every later phase drops excluded
   cells. The excluded share per section is reported, and compared across sections in the young-against-old comparison.
   - Plot: the calls drawn on tissue, to check where each type sits, such as germ cells inside
     tubules and Leydig cells outside them.
   - Plot: a bar chart of cells per call and of how many were excluded.

**Later: reclustering a cell type.** Take the cells of one major type and run the Phase C recipe
on them again, then annotate the subclusters against that type's finer reference. It would
replace step 6 for a type whose fine calls need clusters. Which types have a finer reference:

| Major type | Fine reference | Status |
|---|---|---|
| Germ (every cell called germ) | Green GC1-12 | **Done** by correlation, step 6 |
| Sertoli, SPG subtypes | The collaborator's file, if cleared for use | Open question |
| Leydig, Macrophage, Endothelial, Myoid, InnateLymph, Unknown (Tcf21+) | None yet | Later; other references to find |

## Decisions and open questions

- [ ] Which sections are young and which old, and how many animals in each? (Do the a1, b2, c1,
      d2, 2a, 2b prefixes mark animals?)
- [ ] Green SPG and Sertoli subtypes: the collaborator's file has them, but they are not in the public data. Use them, and if so, where did they come from?
- [ ] Which image channels exist for tubule segmentation?
- [ ] Has the collaborator staged any of these sections? That would be a ready validation set for Phase F.
- [ ] Which ligand-receptor pairs the panel covers
