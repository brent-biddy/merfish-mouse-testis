# MERFISH testis — cell type annotation
Brent Biddy
2026-09-25

- [Setup](#setup)
- [b2r0_cellpose3d](#b2r0_cellpose3d)

# Setup

<details>
<summary>Code</summary>

``` python
import traceback
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import spatialdata
from IPython.display import Markdown, display
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import NullFormatter, ScalarFormatter
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

# The sweep writes twenty; this is the one the document reads.
RESOLUTION = 1.0

# A cluster whose top share falls below this is called Ambiguous, not named on a split vote.
PURITY_THRESHOLD = 0.70

# The template is 16:9 at 10 x 5.625in, minus room for a slide title.
FIGSIZE = (9.0, 4.3)

# Only the columns every sample has: a VPT cell_metadata.csv carries volume alone. n_counts
# is filter_cells' per-cell total, which is what transcript_count would have been.
QC_METRICS = {
    "n_counts": "Transcripts",
    "n_genes": "Genes detected",
    "volume": "Cell volume",
}
# n_genes is capped by the panel, so it spans one decade and a log axis only clutters it.
LOG_METRICS = ("n_counts", "volume")

VIOLIN_INNER = {"box_width": 2, "whis_width": 1, "marker": "o",
                "markersize": 3.0, "markeredgecolor": "white"}

AMBIGUOUS = "Ambiguous"
AMBIGUOUS_COLOR = "#cccccc"

V1_COLUMN = f"leiden_res_{RESOLUTION:.2f}_v1"
V2_COLUMN = f"leiden_res_{RESOLUTION:.2f}_v2"

# The obs columns this document reads, and the step that writes each. corr_<cell type> is
# absent: annotate_celltypes writes those with cell_type_per_cell, named by the reference.
REQUIRED_OBS = {
    "sample": "create_spatialdata or create_spatialdata_cellpose",
    "volume": "create_spatialdata or create_spatialdata_cellpose",
    "n_counts": "cluster_spatialdata_gpu",
    "n_genes": "cluster_spatialdata_gpu",
    V1_COLUMN: "cluster_spatialdata_gpu",
    "cell_type_per_cell": "annotate_celltypes",
}
# A QC metric with no row above would go unchecked.
assert set(QC_METRICS) <= set(REQUIRED_OBS)

# No dpi here: this cell runs after quarto injects the frontmatter's, and would win over it.
plt.rcParams.update({
    "font.size": 8.5,
    "axes.labelsize": 8.5,
    "axes.titlesize": 9.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
})

# The render step stages every samplesheet path flat and names no columns, so globbing is
# how this notebook asks. The glob names a step, not an extension: two steps' outputs stage
# side by side.
ZARRS = sorted(Path(".").glob("*.annotate_celltypes.zarr"))
if not ZARRS:
    raise FileNotFoundError("No *.annotate_celltypes.zarr staged next to the notebook.")

# Keyed on the sample inside each store, not on staging order. A --group_by run names them
# <sample>.<column>.centroids.h5ad, which the same glob takes.
CENTROIDS = {}
for path in sorted(Path(".").glob("*.centroids.h5ad")):
    store = ad.read_h5ad(path)
    CENTROIDS[str(store.obs["sample"].iloc[0])] = store
if not CENTROIDS:
    raise FileNotFoundError("No *.centroids.h5ad staged next to the notebook.")
```

</details>

## Drawing primitives

<details>
<summary>Code</summary>

``` python
def heatmap(ax, frame, cmap="viridis", vmin=None, vmax=None, fmt="{:.2f}"):
    """A labelled matrix, rows and columns taken from the frame's own index."""
    values = frame.to_numpy()
    image = ax.imshow(values, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
    ax.set_xticks(range(frame.shape[1]), frame.columns, rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(frame.shape[0]), frame.index, fontsize=6)

    midpoint = np.nanmean([np.nanmin(values), np.nanmax(values)])
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            ax.text(column, row, fmt.format(values[row, column]),
                    ha="center", va="center", fontsize=4.5,
                    color="black" if values[row, column] > midpoint else "white")
    return image


def log_axis(ax):
    """A log y axis labelled in plain numbers rather than powers of ten."""
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(ScalarFormatter())
    # A sub-two-decade log axis labels its minor ticks 2 x 10^1; the majors are enough.
    ax.yaxis.set_minor_formatter(NullFormatter())


def diagonal_order(frame):
    """Column order that puts each row's strongest column on the diagonal."""
    order = []
    for row in frame.index:
        best = frame.loc[row].idxmax()
        if best not in order:
            order.append(best)
    return order + [column for column in frame.columns if column not in order]


def markdown_table(frame):
    """A pipe table. Hand-rolled: to_markdown needs tabulate, absent from the container."""
    header = "| " + " | ".join(frame.columns) + " |"
    rule = "|" + "|".join("---" for _ in frame.columns) + "|"
    rows = ["| " + " | ".join(str(value) for value in row) + " |" for row in frame.to_numpy()]
    return "\n".join([header, rule] + rows)
```

</details>

## Slides

Each function emits its heading, its caption and its figure, and does no
work on the data beyond what it draws.

<details>
<summary>Code</summary>

``` python
def qc_sample_slide(adata, sample):
    display(Markdown("## QC metrics"))
    display(Markdown(
        "The sample as a whole. Each panel is one per-cell metric drawn as a violin over "
        "every cell, with the box inside marking the quartiles and the median. Transcripts "
        "and cell volume are on log axes; they span orders of magnitude, where genes "
        "detected is capped by the panel. These are the cells that survived the filter in "
        "step 2, not the raw segmentation."
    ))

    frame = adata.obs[list(QC_METRICS)]

    fig, axes = plt.subplots(1, len(QC_METRICS), figsize=(FIGSIZE[0], 3.0),
                             constrained_layout=True)
    for ax, (column, label) in zip(axes, QC_METRICS.items()):
        sns.violinplot(data=frame, y=column, color="#5DB85D", cut=0, density_norm="width",
                       inner="box", inner_kws=VIOLIN_INNER, linewidth=0.5, ax=ax)
        ax.set_ylabel(label, fontsize=9)
        if column in LOG_METRICS:
            log_axis(ax)
        ax.grid(axis="y", color="#b0b0b0", linewidth=0.6)
        ax.set_axisbelow(True)

    fig.suptitle(f"{sample} · {adata.n_obs:,} cells", fontsize=9, x=0.01, ha="left")
    plt.show()
    plt.close(fig)


def cluster_umap_slide(adata, sample):
    display(Markdown("## Clusters on the UMAP"))
    display(Markdown(
        f"Each point is a cell in UMAP space, coloured by its cluster at resolution "
        f"{RESOLUTION:g}. Clusters are numbered by size, so cluster 1 is the largest."
    ))

    order = list(adata.obs[V1_COLUMN].cat.categories)
    # scanpy reads its palette out of uns, keyed on the column it is colouring by.
    adata.uns[f"{V1_COLUMN}_colors"] = sns.color_palette("tab20", len(order)).as_hex()

    fig, ax = plt.subplots(figsize=FIGSIZE)
    sc.pl.umap(adata, color=V1_COLUMN, ax=ax, show=False, frameon=True,
               legend_loc="right margin", title="")
    fig.suptitle(f"{sample} · {adata.n_obs:,} cells · {len(order)} clusters",
                 fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    plt.show()
    plt.close(fig)


def cluster_similarity_slide(similarity, leaf_order, sample):
    display(Markdown("## Cluster similarity"))
    display(Markdown(
        "Spearman correlation between each pair of cluster centroids, over the panel's "
        "genes. The same matrix is drawn twice, both panels in the dendrogram order it "
        "defines. On the left the clusters keep their size-rank numbers, so they are not "
        "in numerical order; on the right the same clusters are renumbered along the "
        "order they now sit in. A position means the same cluster in both, so reading one "
        "against the other gives the mapping."
    ))

    ordered = similarity.loc[leaf_order, leaf_order]
    panels = {
        "Size rank (v1)": list(leaf_order),
        "Dendrogram order (v2)": [str(rank) for rank in range(1, len(leaf_order) + 1)],
    }

    # Not 0: the values sit in a narrow band well above it that a fixed floor would flatten.
    low = float(ordered.to_numpy().min())

    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE, constrained_layout=True)
    for ax, (title, labels) in zip(axes, panels.items()):
        labelled = ordered.copy()
        labelled.index, labelled.columns = labels, labels
        image = heatmap(ax, labelled, vmin=low, vmax=1.0)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Cluster")

    fig.colorbar(image, ax=axes, shrink=0.6, pad=0.02, aspect=25, label="Spearman r")
    fig.suptitle(sample, fontsize=9, x=0.01, ha="left")
    plt.show()
    plt.close(fig)


def qc_violin_slide(adata, sample):
    display(Markdown("## QC metrics by cluster"))
    display(Markdown(
        "Per-cell QC across the clusters, in v2 order. A cluster that correlates cleanly "
        "with a cell type but carries very few transcripts or genes is contamination "
        "rather than that type."
    ))

    order = list(adata.obs[V2_COLUMN].cat.categories)
    frame = adata.obs[list(QC_METRICS) + [V2_COLUMN]]

    fig, axes = plt.subplots(len(QC_METRICS), 1, figsize=(FIGSIZE[0], 5.0), sharex=True)
    for ax, (column, label) in zip(axes, QC_METRICS.items()):
        sns.violinplot(data=frame, x=V2_COLUMN, y=column, order=order, hue=V2_COLUMN,
                       legend=False, cut=0, density_norm="width", linewidth=0.4,
                       palette="tab20", saturation=1, ax=ax)
        ax.set_ylabel(label, fontsize=7)
        ax.set_xlabel("")
        if column in LOG_METRICS:
            log_axis(ax)
    axes[-1].set_xlabel("Cluster (v2)")
    fig.suptitle(sample, fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    plt.show()
    plt.close(fig)


def percell_calls_umap_slide(adata, sample, cell_types, colors):
    display(Markdown("## Per-cell calls on the UMAP"))
    display(Markdown(
        "Every cell coloured by the reference cell type it correlates with most strongly. "
        "This is a call per cell, made without reference to the clustering, so a cluster "
        "of one colour is agreement between two independent views of the data."
    ))

    adata.uns["cell_type_per_cell_colors"] = [colors[t] for t in cell_types]

    fig, ax = plt.subplots(figsize=FIGSIZE)
    sc.pl.umap(adata, color="cell_type_per_cell", ax=ax, show=False, frameon=True,
               legend_loc="right margin", title="")
    fig.suptitle(f"{sample} · {adata.n_obs:,} cells", fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    plt.show()
    plt.close(fig)


def correlation_umap_slide(coords, correlations, sample):
    display(Markdown("## Correlation on the UMAP"))
    display(Markdown(
        "One panel per reference cell type, each cell coloured by how strongly it "
        "correlates with that type. Values are standardized within each cell, so a panel "
        "shows preference rather than how many genes a cell captured."
    ))

    ncols = 4
    nrows = int(np.ceil(correlations.shape[1] / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(FIGSIZE[0], 2.2 * nrows), squeeze=False)
    for ax, cell_type in zip(axes.ravel(), correlations.columns):
        values = correlations[cell_type].to_numpy()
        order = np.argsort(values)          # high-correlating cells drawn on top
        points = ax.scatter(coords[order, 0], coords[order, 1], c=values[order],
                            s=2, cmap="viridis", linewidths=0)
        ax.set_title(cell_type, fontsize=8)
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(points, ax=ax, shrink=0.75, pad=0.02)

    for ax in axes.ravel()[correlations.shape[1]:]:
        fig.delaxes(ax)                     # a hidden axes still reserves its grid cell

    fig.suptitle(sample, fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    plt.show()
    plt.close(fig)


def correlation_heatmap_slide(panels, sample):
    display(Markdown("## Correlation by cluster"))
    display(Markdown(
        "The per-cell correlations averaged over each cluster, drawn three ways: the mean "
        "as it stands, scaled across each row, and scaled down each column. The values are "
        "standardized per cell upstream, so the first panel is a mean z rather than a mean "
        "correlation, and is comparable down a column as well as across a row. Columns are "
        "ordered so each cluster's best match falls on the diagonal.\n\n"
        "**Only the row scaling is a call.** A column's brightest square is the best home "
        "for that cell type whether or not the type is present at all."
    ))

    fig, axes = plt.subplots(1, len(panels), figsize=FIGSIZE, constrained_layout=True)
    for ax, (title, matrix) in zip(axes, panels.items()):
        # No cell annotations: three panels of eleven columns leaves no room for numbers.
        sns.heatmap(matrix, cmap="viridis", linewidths=0.5, linecolor="white",
                    xticklabels=True, yticklabels=True,
                    cbar_kws={"shrink": 0.6, "pad": 0.02}, ax=ax)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("Cell type", fontsize=8)
        ax.set_ylabel("Cluster v1 (v2)", fontsize=8)
        ax.tick_params(axis="x", labelrotation=45, labelsize=6)
        ax.tick_params(axis="y", labelrotation=0, labelsize=6)
        for label in ax.get_xticklabels():
            label.set_ha("right")

    fig.suptitle(sample, fontsize=9, x=0.01, ha="left")
    plt.show()
    plt.close(fig)


def composition_slide(composition, settled, sample):
    display(Markdown("## Call composition by cluster"))
    display(Markdown(
        f"The share of each cluster's cells whose own best match was each cell type, so a "
        f"row sums to one. A boxed square took at least {PURITY_THRESHOLD:.0%} of its "
        f"cluster and settles it; a cluster with no boxed square is starred and called "
        f"{AMBIGUOUS}."
    ))

    labelled = composition.copy()
    labelled.index = [
        f"{v1} ({rank}){'' if settled[v1] else ' *'}"
        for rank, v1 in enumerate(composition.index, start=1)
    ]

    fig, ax = plt.subplots(figsize=FIGSIZE)
    # Three decimals: at two, a share just under the threshold rounds up to it unboxed.
    image = heatmap(ax, labelled, vmin=0.0, vmax=1.0, fmt="{:.3f}")
    for row, cluster in enumerate(composition.index):
        for column, cell_type in enumerate(composition.columns):
            if composition.loc[cluster, cell_type] >= PURITY_THRESHOLD:
                ax.add_patch(Rectangle((column - 0.5, row - 0.5), 1, 1, fill=False,
                                       edgecolor="#d62728", linewidth=1.5))

    for tick, cluster in zip(ax.get_yticklabels(), composition.index):
        if not settled[cluster]:
            tick.set_color("#d62728")
            tick.set_fontweight("bold")

    ax.set_xlabel("Cell type")
    ax.set_ylabel("Cluster v1 (v2)")
    fig.colorbar(image, ax=ax, shrink=0.7, pad=0.02, label="Share of cluster")
    fig.suptitle(sample, fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    plt.show()
    plt.close(fig)


def celltype_umap_slide(adata, sample, order, colors):
    display(Markdown("## Cell type calls on the UMAP"))
    display(Markdown(
        f"The calls the composition heatmap implies. A cluster that concentrated on one "
        f"cell type takes it; one that did not is {AMBIGUOUS}, in grey, rather than named "
        f"on a split vote."
    ))

    adata.uns["cell_type_colors"] = [colors[t] for t in order]

    fig, ax = plt.subplots(figsize=FIGSIZE)
    sc.pl.umap(adata, color="cell_type", ax=ax, show=False, frameon=True,
               legend_loc="right margin", title="")
    fig.suptitle(f"{sample} · {adata.n_obs:,} cells", fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    plt.show()
    plt.close(fig)


def celltype_spatial_slide(boundaries, sample, order, colors):
    display(Markdown("## Cell types on the tissue"))
    display(Markdown(
        "The same calls in tissue coordinates, drawn as the segmented cell boundaries "
        "rather than as points. Testis cell types are radially organised within the "
        "seminiferous tubule, so a correct annotation shows structure here and a wrong one "
        "shows noise. This is the check the embedding cannot give."
    ))

    # Sized to the tissue's own aspect: a fixed one renders it as a sliver in whitespace.
    left, bottom, right, top = boundaries.total_bounds
    aspect = (top - bottom) / (right - left)
    height = min(4.6, 8.0 * aspect)
    fig, ax = plt.subplots(figsize=(height / aspect, height))

    for cell_type in order:
        subset = boundaries[boundaries["cell_type"] == cell_type]
        if subset.empty:
            continue
        subset.plot(ax=ax, color=colors[cell_type], edgecolor="white", linewidth=0.1,
                    label=cell_type)

    ax.set_aspect("equal")
    # micron y increases down the mosaic; matplotlib draws y up, mirroring the tissue.
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])

    # geopandas does not build legend handles from `label`, so they are made here.
    handles = [Patch(facecolor=colors[t], label=t) for t in order
               if (boundaries["cell_type"] == t).any()]
    ax.legend(handles=handles, fontsize=6, frameon=False, loc="center left",
              bbox_to_anchor=(1.0, 0.5))
    fig.suptitle(f"{sample} · {len(boundaries):,} cells", fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    plt.show()
    plt.close(fig)


def flagged_umap_slide(adata, sample, flagged):
    display(Markdown("## Flagged clusters on the UMAP"))
    display(Markdown(
        "The clusters the composition could not settle, each against the rest of the "
        "sample in grey."
    ))

    ncols = min(4, len(flagged))
    nrows = int(np.ceil(len(flagged) / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(FIGSIZE[0], 2.4 * nrows), squeeze=False)
    for ax, cluster in zip(axes.ravel(), flagged):
        # Everything else NaN, which scanpy draws in na_color: grey without a second scatter.
        mask = adata.obs[V1_COLUMN] == cluster
        adata.obs["_flagged"] = pd.Categorical(mask.map({True: cluster}).where(mask))
        adata.uns["_flagged_colors"] = ["#d62728"]

        sc.pl.umap(adata, color="_flagged", ax=ax, show=False, frameon=True,
                   legend_loc=None, na_color="#e6e6e6",
                   title=f"Cluster {cluster} (n={int(mask.sum()):,})")

    for ax in axes.ravel()[len(flagged):]:
        fig.delaxes(ax)

    fig.suptitle(sample, fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    plt.show()
    plt.close(fig)


def resolve_cluster_slide(adata, correlations, sample, cluster, colors):
    display(Markdown(f"## Resolving cluster {cluster}"))
    display(Markdown(
        "On the left the cluster's cells on the embedding, coloured by their own per-cell "
        "call. On the right the same cells placed by their correlation with the two types "
        "most of them call, with the diagonal marking where a cell correlates equally with "
        "both: two arms off it is two populations, one cloud straddling it is one. Along "
        "the bottom, each QC metric against the rest of the sample."
    ))

    in_cluster = (adata.obs[V1_COLUMN] == cluster).to_numpy()
    calls = adata.obs["cell_type_per_cell"]
    present = calls[in_cluster].value_counts()
    present = list(present[present > 0].index)
    contested = present[:2]

    # Two columns per metric, so the umap and scatter above split the row evenly.
    fig, axes = plt.subplot_mosaic(
        [["umap"] * len(QC_METRICS) + ["scatter"] * len(QC_METRICS),
         [metric for metric in QC_METRICS for _ in range(2)]],
        figsize=(FIGSIZE[0], 6.0), height_ratios=[1.5, 1], constrained_layout=True,
    )

    # Everything outside the cluster NaN, drawn in na_color: the sample in grey for free.
    adata.obs["_resolve"] = pd.Categorical(
        calls.astype(str).where(in_cluster), categories=present
    )
    adata.uns["_resolve_colors"] = [colors[t] for t in present]

    # No legend: a margin one breaks the mosaic layout and an on-data one collides. The
    # palette is the slide above's, and the scatter beside names the two types that matter.
    sc.pl.umap(adata, color="_resolve", ax=axes["umap"], show=False, frameon=True,
               legend_loc=None, na_color="#e6e6e6",
               title=f"Cluster {cluster} (n={int(in_cluster.sum()):,})")

    for cell_type in present:
        mask = in_cluster & (calls == cell_type).to_numpy()
        axes["scatter"].scatter(correlations.loc[mask, contested[0]],
                                correlations.loc[mask, contested[1]],
                                s=10, c=colors[cell_type], linewidths=0)
    span = [correlations.loc[in_cluster, contested].to_numpy().min(),
            correlations.loc[in_cluster, contested].to_numpy().max()]
    axes["scatter"].plot(span, span, color="#888888", linewidth=0.8, linestyle="--", zorder=0)
    axes["scatter"].set_title("The two contested cell types", fontsize=8)
    axes["scatter"].set_xlabel(contested[0])
    axes["scatter"].set_ylabel(contested[1])
    sns.despine(ax=axes["scatter"])

    frame = adata.obs[list(QC_METRICS)].copy()
    frame["group"] = np.where(in_cluster, f"Cluster {cluster}", "Other cells")
    order = [f"Cluster {cluster}", "Other cells"]
    for metric, label in QC_METRICS.items():
        ax = axes[metric]
        sns.violinplot(data=frame, x="group", y=metric, order=order, hue="group",
                       hue_order=order, palette=["#d62728", "#bfbfbf"], legend=False,
                       cut=0, density_norm="width", linewidth=0.4, saturation=1, ax=ax)
        ax.set_title(label, fontsize=7)
        ax.set_xlabel("")
        ax.set_ylabel("")
        sns.despine(ax=ax)

    plt.show()
    plt.close(fig)


def annotation_table_slide(summary):
    display(Markdown("## Cluster annotation summary"))
    display(Markdown(
        "One row per cluster: its size, the cell type the composition settled on, the "
        "share that type took, and the runner-up. `margin` is how far the winner sits "
        "above the runner-up. Nothing is filtered on it — it is reported so a thin call "
        "can be spotted."
    ))
    display(Markdown(markdown_table(summary)))
```

</details>

<details>
<summary>Code</summary>

``` python
for path in ZARRS:
    # ── read one sample, with the counts matrix left on disk ──────────────────────────
    # Its own try: a store that will not open has no sample id to head its section with.
    # Uncaught here, it would abort the cell and take every later sample with it.
    try:
        # X stays lazy; no slide reads a value out of it. obs and var cannot: read_lazy
        # returns them as Dataset2D, which strings_to_categoricals and every scanpy
        # plotting call choke on. scverse/anndata#981.
        adata = ad.experimental.read_lazy(f"{path}/tables/table")
        adata.obs = adata.obs.to_memory()
        adata.var = adata.var.to_memory()

        # sc.pl.umap indexes these with a tuple of dimensions, which a lazy array refuses.
        # Only these two: np.asarray flattens obsm's blank-probe frame to nothing.
        for key in ("X_umap", "spatial"):
            adata.obsm[key] = np.asarray(adata.obsm[key])
    except Exception:
        display(Markdown(f"# {path.name}"))
        display(Markdown(f"```\n{traceback.format_exc()}\n```"))
        continue

    missing = [column for column in REQUIRED_OBS if column not in adata.obs]
    if missing:
        display(Markdown(f"# {path.name}"))
        display(Markdown(
            "This store is missing obs columns this report reads:\n\n"
            + "\n".join(f"- `{column}`, written by {REQUIRED_OBS[column]}"
                        for column in missing)
        ))
        del adata
        continue

    # The sample id comes from inside the object, never the staged filename.
    sample = str(adata.obs["sample"].iloc[0])
    display(Markdown(f"# {sample}"))

    try:
        # read_zarr leaves images and points lazy; shapes come into memory, which is what
        # outlines need. More boundaries than cells -- the clustering filter dropped some.
        shapes = spatialdata.read_zarr(path).shapes
        boundaries = shapes[next(key for key in shapes if key.endswith("_polygons"))]

        # ── each cell's correlation with each reference cell type ─────────────────────
        cell_types = list(adata.obs["cell_type_per_cell"].cat.categories)
        correlations = adata.obs[[f"corr_{cell_type}" for cell_type in cell_types]]
        correlations.columns = cell_types

        colors = dict(zip(cell_types, sc.pl.palettes.default_20))
        colors[AMBIGUOUS] = AMBIGUOUS_COLOR   # the absence of a call, not a class of its own

        # ── how alike the clusters are, in gene space ─────────────────────────────────
        # X is summed CP10K, so X / n_cells then log1p is what the reference centroids
        # hold. From the centroids, so the notebook never opens a counts matrix.
        store = CENTROIDS[sample]
        store = store[store.obs["grouping"].astype(str) == V1_COLUMN]
        centroids = pd.DataFrame(
            np.log1p(np.asarray(store.X) / store.obs["n_cells"].to_numpy()[:, None]),
            index=store.obs["group"].astype(str),
            columns=store.var_names,
        )
        similarity = centroids.T.corr(method="spearman")

        # squareform first: linkage picks its input by shape, so a square matrix is read
        # as observations x features and answers a different question with no error.
        # checks=False for the hair off zero floating point leaves on the diagonal;
        # optimal_ordering for the leaf adjacency that makes the block diagonal readable.
        distance = squareform(1 - similarity.to_numpy(), checks=False)
        leaf_order = [similarity.index[leaf] for leaf in
                      leaves_list(linkage(distance, method="average", optimal_ordering=True))]

        # v2 is position down the dendrogram, so neighbours are transcriptionally adjacent.
        # Render-only: it moves with the clustering, which is what v1 is fixed against.
        v2_of = {v1: str(rank) for rank, v1 in enumerate(leaf_order, start=1)}
        adata.obs[V2_COLUMN] = (
            adata.obs[V1_COLUMN].map(v2_of).astype("category")
            .cat.set_categories([str(rank) for rank in range(1, len(leaf_order) + 1)])
        )

        # ── the per-cell correlations averaged within each cluster ────────────────────
        # The mean of annotate_celltypes' per-cell columns, no new correlation. Where the
        # composition below counts argmax votes, this keeps how strongly a cluster leans.
        profiles = correlations.groupby(adata.obs[V1_COLUMN], observed=True).mean()
        profiles = profiles.loc[leaf_order]
        profiles.index = [f"{v1} ({rank})" for rank, v1 in enumerate(leaf_order, start=1)]

        by_cluster = profiles.sub(profiles.min(axis=1), axis=0).div(
            profiles.max(axis=1) - profiles.min(axis=1), axis=0)
        by_type = profiles.sub(profiles.min(axis=0), axis=1).div(
            profiles.max(axis=0) - profiles.min(axis=0), axis=1)

        type_order = diagonal_order(by_cluster)
        correlation_panels = {
            "Mean per-cell z": profiles[type_order],
            "Scaled per cluster": by_cluster[type_order],
            "Scaled per cell type": by_type[type_order],
        }

        # ── what the per-cell calls compose to inside each cluster ────────────────────
        counts = pd.crosstab(adata.obs[V1_COLUMN], adata.obs["cell_type_per_cell"])
        composition = counts.div(counts.sum(axis=1), axis=0).loc[leaf_order]
        # Its own ordering: a cluster's largest share and strongest mean correlation differ.
        composition = composition[diagonal_order(composition)]

        settled = composition.max(axis=1) >= PURITY_THRESHOLD
        calls = {
            cluster: (composition.loc[cluster].idxmax() if is_settled else AMBIGUOUS)
            for cluster, is_settled in settled.items()
        }
        flagged = [cluster for cluster in leaf_order if not settled[cluster]]

        called_types = sorted({call for call in calls.values() if call != AMBIGUOUS})
        cell_type_order = called_types + [AMBIGUOUS]          # Ambiguous last in the legend
        adata.obs["cell_type"] = (
            adata.obs[V1_COLUMN].map(calls).astype("category")
            .cat.set_categories(cell_type_order)
        )

        boundaries = boundaries.join(adata.obs["cell_type"], how="inner")

        summary_rows = []
        for v1 in leaf_order:
            shares = composition.loc[v1].sort_values(ascending=False)
            summary_rows.append({
                "cluster_v1": v1,
                "cluster_v2": v2_of[v1],
                "n_cells": f"{int(counts.loc[v1].sum()):,}",
                "call": calls[v1],
                "top1": shares.index[0],
                "share": f"{shares.iloc[0]:.3f}",
                "top2": shares.index[1],
                "margin": f"{shares.iloc[0] - shares.iloc[1]:.3f}",
                "flagged": "*" if v1 in flagged else "",
            })
        summary = pd.DataFrame(summary_rows)
    except Exception:
        display(Markdown(f"```\n{traceback.format_exc()}\n```"))
        del adata
        continue

    # ── the slides, in the order they appear ──────────────────────────────────────────
    slides = [
        (qc_sample_slide, (adata, sample)),
        (cluster_umap_slide, (adata, sample)),
        (cluster_similarity_slide, (similarity, leaf_order, sample)),
        (qc_violin_slide, (adata, sample)),
        (percell_calls_umap_slide, (adata, sample, cell_types, colors)),
        (correlation_umap_slide, (adata.obsm["X_umap"], correlations, sample)),
        (correlation_heatmap_slide, (correlation_panels, sample)),
        (composition_slide, (composition, settled, sample)),
        (celltype_umap_slide, (adata, sample, cell_type_order, colors)),
        (celltype_spatial_slide, (boundaries, sample, cell_type_order, colors)),
    ]
    if flagged:
        slides.append((flagged_umap_slide, (adata, sample, flagged)))
        slides += [(resolve_cluster_slide, (adata, correlations, sample, cluster, colors))
                   for cluster in flagged]
    slides.append((annotation_table_slide, (summary,)))

    # Per slide: each draws its own figure, and one failing says nothing about the others.
    for slide, arguments in slides:
        try:
            slide(*arguments)
        except Exception:
            plt.close()
            display(Markdown(f"```\n{traceback.format_exc()}\n```"))

    del adata
```

</details>

# b2r0_cellpose3d

## QC metrics

The sample as a whole. Each panel is one per-cell metric drawn as a
violin over every cell, with the box inside marking the quartiles and
the median. Transcripts and cell volume are on log axes; they span
orders of magnitude, where genes detected is capped by the panel. These
are the cells that survived the filter in step 2, not the raw
segmentation.

![](celltype_report_files/figure-commonmark/cell-5-output-4.png)

## Clusters on the UMAP

Each point is a cell in UMAP space, coloured by its cluster at
resolution 1. Clusters are numbered by size, so cluster 1 is the
largest.

![](celltype_report_files/figure-commonmark/cell-5-output-7.png)

## Cluster similarity

Spearman correlation between each pair of cluster centroids, over the
panel’s genes. The same matrix is drawn twice, both panels in the
dendrogram order it defines. On the left the clusters keep their
size-rank numbers, so they are not in numerical order; on the right the
same clusters are renumbered along the order they now sit in. A position
means the same cluster in both, so reading one against the other gives
the mapping.

![](celltype_report_files/figure-commonmark/cell-5-output-10.png)

## QC metrics by cluster

Per-cell QC across the clusters, in v2 order. A cluster that correlates
cleanly with a cell type but carries very few transcripts or genes is
contamination rather than that type.

![](celltype_report_files/figure-commonmark/cell-5-output-13.png)

## Per-cell calls on the UMAP

Every cell coloured by the reference cell type it correlates with most
strongly. This is a call per cell, made without reference to the
clustering, so a cluster of one colour is agreement between two
independent views of the data.

![](celltype_report_files/figure-commonmark/cell-5-output-16.png)

## Correlation on the UMAP

One panel per reference cell type, each cell coloured by how strongly it
correlates with that type. Values are standardized within each cell, so
a panel shows preference rather than how many genes a cell captured.

![](celltype_report_files/figure-commonmark/cell-5-output-19.png)

## Correlation by cluster

The per-cell correlations averaged over each cluster, drawn three ways:
the mean as it stands, scaled across each row, and scaled down each
column. The values are standardized per cell upstream, so the first
panel is a mean z rather than a mean correlation, and is comparable down
a column as well as across a row. Columns are ordered so each cluster’s
best match falls on the diagonal.

**Only the row scaling is a call.** A column’s brightest square is the
best home for that cell type whether or not the type is present at all.

![](celltype_report_files/figure-commonmark/cell-5-output-22.png)

## Call composition by cluster

The share of each cluster’s cells whose own best match was each cell
type, so a row sums to one. A boxed square took at least 70% of its
cluster and settles it; a cluster with no boxed square is starred and
called Ambiguous.

![](celltype_report_files/figure-commonmark/cell-5-output-25.png)

## Cell type calls on the UMAP

The calls the composition heatmap implies. A cluster that concentrated
on one cell type takes it; one that did not is Ambiguous, in grey,
rather than named on a split vote.

![](celltype_report_files/figure-commonmark/cell-5-output-28.png)

## Cell types on the tissue

The same calls in tissue coordinates, drawn as the segmented cell
boundaries rather than as points. Testis cell types are radially
organised within the seminiferous tubule, so a correct annotation shows
structure here and a wrong one shows noise. This is the check the
embedding cannot give.

![](celltype_report_files/figure-commonmark/cell-5-output-31.png)

## Flagged clusters on the UMAP

The clusters the composition could not settle, each against the rest of
the sample in grey.

![](celltype_report_files/figure-commonmark/cell-5-output-34.png)

## Resolving cluster 16

On the left the cluster’s cells on the embedding, coloured by their own
per-cell call. On the right the same cells placed by their correlation
with the two types most of them call, with the diagonal marking where a
cell correlates equally with both: two arms off it is two populations,
one cloud straddling it is one. Along the bottom, each QC metric against
the rest of the sample.

![](celltype_report_files/figure-commonmark/cell-5-output-37.png)

## Resolving cluster 14

On the left the cluster’s cells on the embedding, coloured by their own
per-cell call. On the right the same cells placed by their correlation
with the two types most of them call, with the diagonal marking where a
cell correlates equally with both: two arms off it is two populations,
one cloud straddling it is one. Along the bottom, each QC metric against
the rest of the sample.

![](celltype_report_files/figure-commonmark/cell-5-output-40.png)

## Resolving cluster 18

On the left the cluster’s cells on the embedding, coloured by their own
per-cell call. On the right the same cells placed by their correlation
with the two types most of them call, with the diagonal marking where a
cell correlates equally with both: two arms off it is two populations,
one cloud straddling it is one. Along the bottom, each QC metric against
the rest of the sample.

![](celltype_report_files/figure-commonmark/cell-5-output-43.png)

## Resolving cluster 11

On the left the cluster’s cells on the embedding, coloured by their own
per-cell call. On the right the same cells placed by their correlation
with the two types most of them call, with the diagonal marking where a
cell correlates equally with both: two arms off it is two populations,
one cloud straddling it is one. Along the bottom, each QC metric against
the rest of the sample.

![](celltype_report_files/figure-commonmark/cell-5-output-46.png)

## Resolving cluster 13

On the left the cluster’s cells on the embedding, coloured by their own
per-cell call. On the right the same cells placed by their correlation
with the two types most of them call, with the diagonal marking where a
cell correlates equally with both: two arms off it is two populations,
one cloud straddling it is one. Along the bottom, each QC metric against
the rest of the sample.

![](celltype_report_files/figure-commonmark/cell-5-output-49.png)

## Resolving cluster 15

On the left the cluster’s cells on the embedding, coloured by their own
per-cell call. On the right the same cells placed by their correlation
with the two types most of them call, with the diagonal marking where a
cell correlates equally with both: two arms off it is two populations,
one cloud straddling it is one. Along the bottom, each QC metric against
the rest of the sample.

![](celltype_report_files/figure-commonmark/cell-5-output-52.png)

## Resolving cluster 19

On the left the cluster’s cells on the embedding, coloured by their own
per-cell call. On the right the same cells placed by their correlation
with the two types most of them call, with the diagonal marking where a
cell correlates equally with both: two arms off it is two populations,
one cloud straddling it is one. Along the bottom, each QC metric against
the rest of the sample.

![](celltype_report_files/figure-commonmark/cell-5-output-55.png)

## Resolving cluster 6

On the left the cluster’s cells on the embedding, coloured by their own
per-cell call. On the right the same cells placed by their correlation
with the two types most of them call, with the diagonal marking where a
cell correlates equally with both: two arms off it is two populations,
one cloud straddling it is one. Along the bottom, each QC metric against
the rest of the sample.

![](celltype_report_files/figure-commonmark/cell-5-output-58.png)

## Resolving cluster 12

On the left the cluster’s cells on the embedding, coloured by their own
per-cell call. On the right the same cells placed by their correlation
with the two types most of them call, with the diagonal marking where a
cell correlates equally with both: two arms off it is two populations,
one cloud straddling it is one. Along the bottom, each QC metric against
the rest of the sample.

![](celltype_report_files/figure-commonmark/cell-5-output-61.png)

## Resolving cluster 1

On the left the cluster’s cells on the embedding, coloured by their own
per-cell call. On the right the same cells placed by their correlation
with the two types most of them call, with the diagonal marking where a
cell correlates equally with both: two arms off it is two populations,
one cloud straddling it is one. Along the bottom, each QC metric against
the rest of the sample.

![](celltype_report_files/figure-commonmark/cell-5-output-64.png)

## Cluster annotation summary

One row per cluster: its size, the cell type the composition settled on,
the share that type took, and the runner-up. `margin` is how far the
winner sits above the runner-up. Nothing is filtered on it — it is
reported so a thin call can be spotted.

| cluster_v1 | cluster_v2 | n_cells | call | top1 | share | top2 | margin | flagged |
|----|----|----|----|----|----|----|----|----|
| 16 | 1 | 4,544 | Ambiguous | Myoid | 0.478 | ImmLeydig | 0.236 | \* |
| 14 | 2 | 5,223 | Ambiguous | ImmLeydig | 0.351 | f-Pericyte | 0.108 | \* |
| 18 | 3 | 2,312 | Ambiguous | Endothelial | 0.443 | f-Pericyte | 0.176 | \* |
| 17 | 4 | 3,794 | Spermatogonia | Spermatogonia | 0.709 | Spermatocyte | 0.593 |  |
| 11 | 5 | 6,721 | Ambiguous | Spermatogonia | 0.485 | Spermatocyte | 0.261 | \* |
| 13 | 6 | 5,933 | Ambiguous | ImmLeydig | 0.260 | f-Pericyte | 0.043 | \* |
| 8 | 7 | 7,647 | Spermatocyte | Spermatocyte | 0.726 | Spermatogonia | 0.527 |  |
| 5 | 8 | 11,478 | Spermatocyte | Spermatocyte | 0.789 | Elongating | 0.714 |  |
| 3 | 9 | 12,247 | Spermatocyte | Spermatocyte | 0.780 | Elongating | 0.636 |  |
| 10 | 10 | 6,951 | Spermatocyte | Spermatocyte | 0.866 | Elongating | 0.771 |  |
| 15 | 11 | 4,962 | Ambiguous | RoundSpermatid | 0.461 | Elongating | 0.112 | \* |
| 19 | 12 | 1,071 | Ambiguous | Spermatocyte | 0.531 | RoundSpermatid | 0.201 | \* |
| 7 | 13 | 8,184 | Spermatocyte | Spermatocyte | 0.955 | RoundSpermatid | 0.920 |  |
| 9 | 14 | 7,306 | Spermatocyte | Spermatocyte | 0.707 | RoundSpermatid | 0.430 |  |
| 2 | 15 | 13,044 | RoundSpermatid | RoundSpermatid | 0.784 | Spermatocyte | 0.660 |  |
| 6 | 16 | 10,580 | Ambiguous | Elongating | 0.527 | RoundSpermatid | 0.059 | \* |
| 12 | 17 | 6,378 | Ambiguous | Elongating | 0.632 | RoundSpermatid | 0.266 | \* |
| 1 | 18 | 17,308 | Ambiguous | Elongating | 0.567 | RoundSpermatid | 0.134 | \* |
| 4 | 19 | 11,685 | Elongating | Elongating | 0.713 | RoundSpermatid | 0.507 |  |
