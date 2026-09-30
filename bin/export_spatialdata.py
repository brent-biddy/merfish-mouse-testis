#!/usr/bin/env python3
"""
export_spatialdata.py - Put a store and a table back together as one SpatialData zarr.

Every step after create_spatialdata writes only the table, leaving the images, transcripts
and boundaries in the store create_spatialdata wrote. This writes the two as one
self-contained store, for sharing or for a tool that opens a SpatialData directly. It is
the only step after create_spatialdata that copies the images, so run it when a single
object is wanted rather than as part of every analysis.

The table names the shapes element it annotates. SpatialData only warns when an object is
built or read with a table whose element is missing, and says nothing when the table is
assigned to an object already open, as here -- so that is checked, as is that the table's
sample is the one asked for.

Writes <outdir>/<sample>.export_spatialdata.zarr plus a timing TSV.

Usage:
    export_spatialdata.py --sample testis_01 \\
        --zarr_store results/testis_01/create_spatialdata/testis_01.create_spatialdata.zarr \\
        --table_path results/testis_01/annotate_celltypes/testis_01.annotate_celltypes.h5ad \\
        --outdir results/testis_01/export_spatialdata
"""

import argparse
from pathlib import Path

import anndata as ad
import spatialdata

from timer import timer, timing_summary


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine a SpatialData store and its table into one zarr"
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="Sample identifier",
    )
    parser.add_argument(
        "--zarr_store",
        required=True,
        help="Store from create_spatialdata or create_spatialdata_cellpose",
    )
    parser.add_argument(
        "--table_path",
        required=True,
        help="Table .h5ad from any step, e.g. annotate_celltypes",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory to write <sample>.export_spatialdata.zarr into (default: current directory)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_path = outdir / f"{args.sample}.export_spatialdata.zarr"

    print(f"Sample:  {args.sample}")
    print(f"Store:   {args.zarr_store}")
    print(f"Table:   {args.table_path}")
    print(f"Output:  {output_path}")

    with timer("Read store"):
        sdata = spatialdata.read_zarr(args.zarr_store)

    with timer("Read table"):
        table = ad.read_h5ad(args.table_path)

    samples = sorted(table.obs["sample"].astype(str).unique())
    if samples != [args.sample]:
        raise ValueError(f"{args.table_path} holds sample(s) {samples}, not {args.sample}.")

    annotated = table.uns["spatialdata_attrs"]["region"]
    annotated = [annotated] if isinstance(annotated, str) else list(annotated)
    missing = [name for name in annotated if name not in sdata.shapes]
    if missing:
        raise ValueError(
            f"The table annotates {missing}, which {args.zarr_store} does not hold; it has "
            f"shapes {sorted(sdata.shapes)}. The store and the table are not from one sample."
        )

    if "table" in sdata.tables:
        print("The store already held a table; replacing it with the one given.")
    sdata.tables["table"] = table

    print(f"Table:   {table.n_obs:,} cells x {table.n_vars:,} genes, annotating {annotated}")

    with timer("Write zarr"):
        sdata.write(output_path, overwrite=True)

    timing_summary(outdir / f"{args.sample}.export_spatialdata.timing.tsv")


if __name__ == "__main__":
    main()
