#!/usr/bin/env python3
"""
cluster_spatialdata_gpu.py - GPU-accelerated QC, normalize, and cluster a SpatialData zarr.

Uses rapids-singlecell for the compute-heavy steps (filtering, normalization, PCA,
neighbors, UMAP, and a sweep of Leiden clusterings). Data is moved back to CPU before
zarr I/O.

There is no highly-variable-gene selection: a MERFISH panel is a few hundred curated
markers, so every gene is used.

Each swept resolution leaves two obs columns, e.g. leiden_res_0.10_v0 and leiden_res_0.10_v1.
v1 is the size ranking and is what downstream steps mean by a cluster id.

Writes <outdir>/<sample>.cluster_spatialdata_gpu.zarr plus a timing TSV.

Usage:
    cluster_spatialdata_gpu.py --sample testis_01 \\
        --path results/testis_01/create_spatialdata/testis_01.zarr \\
        --outdir results/testis_01/cluster_spatialdata_gpu
"""

import argparse
from pathlib import Path

import rapids_singlecell as rsc
import spatialdata

from timer import timer, timing_summary

# Leiden resolutions to sweep.
RESOLUTIONS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
               1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]

SCALE_MAX_VALUE = 10
N_PCS = 30
NEIGHBORS_METRIC = "cosine"


def relabel_by_size(labels):
    """Renumber Leiden labels 1..k by descending cluster size, so cluster "1" is the largest.

    This is the v0 -> v1 step: `labels` is leiden's own output and the return value is its
    size-ranked sibling. Both are kept, so a v1 id can always be traced back.

    Returns an ordered categorical, so plots and tables sort 2 before 10.
    """
    largest_first = labels.value_counts().index      # value_counts sorts descending by size
    ranked_categories = [str(rank) for rank in range(1, len(largest_first) + 1)]
    size_rank = dict(zip(largest_first, ranked_categories))

    return (
        labels.map(size_rank)
        .astype("category")
        .cat.set_categories(ranked_categories, ordered=True)
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="GPU-accelerated clustering of a SpatialData zarr"
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="Sample identifier",
    )
    parser.add_argument(
        "--path",
        required=True,
        help="Path to input SpatialData zarr",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory to write <sample>.cluster_spatialdata_gpu.zarr into (default: current directory)",
    )
    parser.add_argument(
        "--min_counts",
        type=int,
        default=20,
        help="Drop cells with fewer than this many transcripts (default: %(default)s).",
    )
    parser.add_argument(
        "--max_counts_quantile",
        type=float,
        default=0,
        help="Drop cells above this quantile of transcript_count, to remove doublets and "
        "segmentation merges. 0 disables the cut (default: %(default)s).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_path = outdir / f"{args.sample}.cluster_spatialdata_gpu.zarr"

    print(f"Sample:  {args.sample}")
    print(f"Input:   {args.path}")
    print(f"Output:  {output_path}")
    print(f"Res:     {', '.join(f'{r:g}' for r in RESOLUTIONS)}")
    print(f"Filter:  min_counts={args.min_counts}, "
          f"max_counts_quantile={args.max_counts_quantile:g}")

    with timer("Read zarr"):
        sdata = spatialdata.read_zarr(args.path)

    with timer("Extract table"):
        adata = sdata.tables["table"].copy()

    print(f"Table:   {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Record the parameters this step applied, so downstream code reads them from the
    # object rather than restating them.
    adata.uns["cluster_min_counts"] = args.min_counts
    adata.uns["cluster_max_counts_quantile"] = args.max_counts_quantile
    adata.uns["cluster_resolutions"] = RESOLUTIONS

    # merscope() returns int64 counts, and the rapids QC kernels return zeros for int64
    # input instead of erroring — every cell then reads as empty and gets filtered out.
    adata.X = adata.X.astype("float32")

    with timer("Move to GPU"):
        rsc.get.anndata_to_GPU(adata)

    with timer("Filter"):
        n_before = adata.n_obs
        # Quantile on the unfiltered table, before min_counts removes the low tail.
        max_counts = (
            float(adata.obs["transcript_count"].quantile(args.max_counts_quantile))
            if args.max_counts_quantile else None
        )
        rsc.pp.filter_cells(adata, min_counts=args.min_counts)
        if max_counts is not None:
            rsc.pp.filter_cells(adata, max_counts=max_counts)
            print(f"Upper cut at q{args.max_counts_quantile:g} = {max_counts:,.0f} transcripts.")

    print(f"Filtered {n_before - adata.n_obs:,} low-quality cells.")
    print(f"Retained {adata.n_obs:,} cells x {adata.n_vars:,} genes.")

    with timer("Normalize"):
        adata.layers["counts"] = adata.X.copy()
        # No target_sum, so the default applies: the median pre-normalization cell total.
        # That makes lognorm a per-sample scale -- compare across samples from counts.
        rsc.pp.normalize_total(adata, inplace=True)
        rsc.pp.log1p(adata)
        # Preserved before scaling overwrites X; downstream annotation reads this.
        adata.layers["lognorm"] = adata.X.copy()

    with timer("Scale"):
        rsc.pp.scale(adata, zero_center=False, max_value=SCALE_MAX_VALUE)

    with timer("PCA"):
        rsc.pp.pca(adata, n_comps=N_PCS, random_state=0)

    with timer("Neighbors"):
        rsc.pp.neighbors(adata, metric=NEIGHBORS_METRIC, random_state=0)

    with timer("UMAP"):
        rsc.tl.umap(adata, random_state=0)

    # Sweep resolutions rather than committing to one.
    for res in RESOLUTIONS:
        leiden_key = f"leiden_res_{res:.2f}_v0"
        ranked_key = f"leiden_res_{res:.2f}_v1"
        with timer(f"Leiden res={res:g}"):
            rsc.tl.leiden(
                adata,
                resolution=res,
                key_added=leiden_key,
                random_state=0,
            )
            adata.obs[ranked_key] = relabel_by_size(adata.obs[leiden_key])

    # rapids-singlecell keeps arrays on GPU; zarr I/O needs them back on the host.
    with timer("Move to CPU"):
        rsc.get.anndata_to_CPU(adata)

    with timer("Write zarr"):
        sdata.tables["table"] = adata
        sdata.write(output_path, overwrite=True)

    print(f"Written to {output_path}")

    timing_summary(outdir / f"{args.sample}.cluster_spatialdata_gpu.timing.tsv")


if __name__ == "__main__":
    main()
