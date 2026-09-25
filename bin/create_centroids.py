#!/usr/bin/env python3
"""
create_centroids.py - Build per-cluster centroids from a clustered SpatialData zarr.

Sums both expression layers over an obs column, writing one row per group of cells
instead of one per cell:

    X                 each cell's CP10K profile, summed over the group
    layers["counts"]  raw counts, summed over the group

Stored as sums, not means, so clusters pool: any union's profile is the row-wise sum of
its members, and n_cells sums with it. The reference centroids in assets/reference are
ln(mean + 1), so the comparable value here is log1p(X / n_cells).

Groups by every _v1 column by default, or by one named obs column with --group_by.
Requires layers["counts"], which cluster_spatialdata_gpu.py writes.

Writes <outdir>/<sample>.centroids.h5ad plus a timing TSV, or
<outdir>/<sample>.<column>.centroids.* for a --group_by run, so the two can sit side
by side.

Usage:
    create_centroids.py --sample testis_01 \\
        --path results/testis_01/cluster_spatialdata_gpu/testis_01.zarr \\
        --outdir results/testis_01/create_centroids
"""

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
import spatialdata

from timer import timer, timing_summary


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build per-cluster centroids from a clustered SpatialData zarr"
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="Sample identifier",
    )
    parser.add_argument(
        "--path",
        required=True,
        help="Clustered SpatialData zarr from cluster_spatialdata_gpu",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory to write <sample>.centroids.h5ad into (default: current directory)",
    )
    parser.add_argument(
        "--group_by",
        default=None,
        help="obs column to sum over. Omit to sum over every _v1 column in the sweep.",
    )
    return parser.parse_args()


def get_centroids(adata, column):
    """Sum both layers over one obs column, returning a groups x genes AnnData."""
    # sc.get.aggregate takes one layer at a time, so two calls.
    cp10k = sc.get.aggregate(adata, by=column, func="sum", layer="cp10k")
    counts = sc.get.aggregate(adata, by=column, func="sum", layer="counts")

    out = ad.AnnData(
        X=cp10k.layers["sum"],
        obs=pd.DataFrame(index=cp10k.obs_names),
        var=pd.DataFrame(index=cp10k.var_names),
        layers={"counts": counts.layers["sum"]},
    )

    out.obs["grouping"] = column
    out.obs["group"] = out.obs_names
    out.obs["n_cells"] = cp10k.obs["n_obs_aggregated"].to_numpy()
    out.obs_names = [f"{column}_{group}" for group in out.obs["group"]]

    return out


def get_grouping_columns(adata, group_by):
    """The obs columns to build centroids for, validated before the expensive CP10K pass."""
    if group_by is not None:
        if group_by not in adata.obs:
            raise ValueError(
                f"--group_by {group_by!r} is not an obs column of this zarr. This step "
                f"groups by a column an earlier step wrote; it does not compute labels. "
                f"Found: {sorted(adata.obs.columns)}"
            )
        return [group_by]

    # _v1 only: the _v0 siblings would cover every resolution twice.
    # Already in sweep order, so no sort.
    columns = [column for column in adata.obs if column.endswith("_v1")]
    if not columns:
        raise ValueError(
            f"no _v1 columns in obs — was this zarr written by cluster_spatialdata_gpu? "
            f"Found: {sorted(adata.obs.columns)}"
        )
    return columns


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # A --group_by run publishes beside a sweep run, so the column goes in the
    # name: <sample>.<column>.centroids, not <sample>.centroids.<column>, so one
    # glob takes both.
    stem = f"{args.sample}.{args.group_by}.centroids" if args.group_by else f"{args.sample}.centroids"
    output_path = outdir / f"{stem}.h5ad"

    print(f"Sample:  {args.sample}")
    print(f"Input:   {args.path}")
    print(f"Output:  {output_path}")
    print(f"Group:   {args.group_by or 'leiden sweep'}")

    with timer("Read zarr"):
        sdata = spatialdata.read_zarr(args.path)

    with timer("Extract table"):
        adata = sdata.tables["table"]

    print(f"Table:   {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    if "counts" not in adata.layers:
        raise KeyError(
            f"{args.sample}: the table has no layers['counts'], so centroids cannot be "
            f"built. cluster_spatialdata_gpu.py writes it; re-run that step over this "
            f"sample. Layers present: {sorted(adata.layers)}"
        )

    columns = get_grouping_columns(adata, args.group_by)
    print(f"Centroids for {len(columns)} column(s): {', '.join(columns)}")

    with timer("CP10K"):
        # From layers["counts"], not expm1(X): scaling overwrote X, so it no longer
        # holds expression.
        adata.layers["cp10k"] = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, layer="cp10k", target_sum=1e4)

    centroid_list = []
    for column in columns:
        with timer(f"Centroids {column}"):
            centroid_list.append(get_centroids(adata, column))

    with timer("Assemble"):
        centroids = ad.concat(centroid_list, axis=0, join="outer")
        centroids.obs["sample"] = args.sample

        # anndata leaves a string column as object only when every value is unique, so
        # `group` varies: categorical from a sweep, object from a --group_by run.
        for column in ("grouping", "group", "sample"):
            centroids.obs[column] = pd.Categorical(centroids.obs[column].astype(str))

    print(f"Assembled {centroids.n_obs:,} group rows x {centroids.n_vars:,} genes.")

    with timer("Write h5ad"):
        centroids.write_h5ad(output_path)

    timing_summary(outdir / f"{stem}.timing.tsv")


if __name__ == "__main__":
    main()
