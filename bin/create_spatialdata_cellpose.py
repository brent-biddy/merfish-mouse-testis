#!/usr/bin/env python3
"""
create_spatialdata_cellpose.py - Convert a cellpose re-segmentation of a MERSCOPE region
to a SpatialData Zarr store.

Reads two directories. The MERSCOPE region directory supplies the mosaic images, which the
re-segmentation does not touch. The VPT output directory supplies everything it does: the
count matrix, the cell metadata, the boundary polygons, and the detected transcripts, which
VPT rewrites with a cell_id column naming the cell each one falls in. The transcripts are
staged in over the region's own copy, which has no such column.

Both are named on the command line. Which segmentation a store came from is not
recoverable from its contents, so the VPT directory is recorded into the table.

The store is the same shape as the one create_spatialdata.py writes, so the clustering,
annotation and report steps read it unchanged.

Writes <outdir>/<sample>.create_spatialdata_cellpose.zarr plus a timing TSV.

Usage:
    create_spatialdata_cellpose.py --sample b2r0_cellpose \\
        --path data/raw/MsTestis_B2/region_0 \\
        --vpt_path data/cellpose/MsTestis_B2/region_0/Cellpose \\
        --outdir results/b2r0_cellpose/create_spatialdata_cellpose
"""

import argparse
from pathlib import Path

import pandas as pd
from spatialdata_io import merscope

from timer import timer, timing_summary

# merscope()'s z_layers default, pinned so a library change cannot move the plane silently.
# It selects the image alone: one per-stain mosaic TIFF set, read into the element
# <slide>_<region>_z<n>. Recorded into the table below as provenance; nothing reads it back.
Z_LAYER = 3

# Keyed by the names merscope()'s vpt_outputs dict takes. These are the Vizgen delivery's
# names, which prep_cellpose_vpt.py writes too. VPT itself prefixes only the boundaries --
# a stock run writes cell_by_gene.csv and cell_metadata.csv unprefixed, and would have to
# be named here. The dict form is needed regardless: the reordering below hands merscope()
# a rewritten metadata path, which the directory form has no way to take.
VPT_FILES = {
    "cell_by_gene": "cellpose_cell_by_gene.csv",
    "cell_metadata": "cellpose_cell_metadata.csv",
    "cell_boundaries": "cellpose_micron_space.parquet",
}

# VPT rewrites the transcripts too, unprefixed, adding a cell_id column the region's own
# copy does not have. Not in VPT_FILES: merscope()'s dict takes only the three above.
TRANSCRIPTS_FILE = "detected_transcripts.csv"


def vpt_outputs(vpt_dir, staging_dir):
    """Return the vpt_outputs dict merscope() reads, reordering the metadata if needed.

    merscope() hands the counts and the metadata straight to AnnData, which requires their
    indexes to match exactly; VPT does not guarantee it writes them in the same order. A
    reordered copy is written to staging_dir rather than over the input.

    A missing file is an error here, whichever of the three it is, because merscope() only
    warns. Without the boundaries it loads no polygons, leaving a store that clusters and
    annotates perfectly well and has nothing to draw a tissue figure from; without either
    CSV it loads no table at all. Neither is a thing to discover several steps later.
    """
    paths = {key: vpt_dir / name for key, name in VPT_FILES.items()}

    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(
            f"Missing VPT output(s) in {vpt_dir}: {', '.join(missing)}. Expected the names a "
            f"Vizgen cellpose delivery uses; a stock VPT run writes the two CSVs unprefixed, "
            f"and a watershed run names its boundaries watershed_micron_space.parquet."
        )

    # Only the ids are needed, and cellpose_cell_by_gene.csv is hundreds of MB of counts.
    # dtype={0: str} rather than dtype=str: a scalar dtype does not reach index_col.
    counts_index = pd.read_csv(
        paths["cell_by_gene"], index_col=0, usecols=[0], dtype={0: str}
    ).index
    metadata = pd.read_csv(paths["cell_metadata"], index_col=0, dtype={0: str})

    if not counts_index.equals(metadata.index):
        print(f"Reordering cell metadata to match the counts ({len(counts_index):,} cells).")
        staging_dir.mkdir(parents=True, exist_ok=True)
        paths["cell_metadata"] = staging_dir / VPT_FILES["cell_metadata"]
        metadata.loc[counts_index].to_csv(paths["cell_metadata"])

    return paths


def staged_region_dir(region_dir, vpt_dir, staging_dir):
    """Return a region directory merscope() can read, staging one with VPT's transcripts.

    One thing it fixes. The transcripts: merscope() takes them from the region directory
    and vpt_outputs overrides only the counts, metadata and boundaries, so VPT's own copy —
    the same rows plus a cell_id column — is reachable only by staging it in over the
    region's. Always stages, unlike create_spatialdata.py's, since that copy is always the
    one wanted.
    """
    vpt_transcripts = vpt_dir / TRANSCRIPTS_FILE
    if not vpt_transcripts.exists():
        raise FileNotFoundError(
            f"No {TRANSCRIPTS_FILE} in {vpt_dir}. merscope() takes transcripts from the "
            f"region directory, whose copy has no cell_id, so the store's points would "
            f"carry no cell assignment while its cells are cellpose cells."
        )

    # Symlink the rest of the region — the mosaic images are far too large to copy — and
    # link VPT's transcripts in over the region's own, so neither input is modified.
    staging_dir.mkdir(parents=True, exist_ok=True)
    for entry in region_dir.iterdir():
        if entry.name != TRANSCRIPTS_FILE:
            # Replaced rather than skipped: exists() follows the link, so one left by an
            # earlier run pointing somewhere gone reads as absent and symlink_to then fails.
            link = staging_dir / entry.name
            link.unlink(missing_ok=True)
            link.symlink_to(entry.resolve())

    link = staging_dir / TRANSCRIPTS_FILE
    link.unlink(missing_ok=True)
    link.symlink_to(vpt_transcripts.resolve())

    return staging_dir


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a cellpose re-segmentation of a MERSCOPE region to a "
                    "SpatialData Zarr store"
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="Sample identifier",
    )
    parser.add_argument(
        "--path",
        required=True,
        help="MERSCOPE region output directory, for the images and transcripts",
    )
    parser.add_argument(
        "--vpt_path",
        required=True,
        help="VPT cellpose output directory, for the cells, counts and boundaries",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory to write <sample>.create_spatialdata_cellpose.zarr into (default: current directory)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_path = outdir / f"{args.sample}.create_spatialdata_cellpose.zarr"

    print(f"Sample:  {args.sample}")
    print(f"Input:   {args.path}")
    print(f"VPT:     {args.vpt_path}")
    print(f"Output:  {output_path}")
    print(f"Z-layer: {Z_LAYER}")

    raw_region_dir = Path(args.path)
    vpt_dir = Path(args.vpt_path)
    paths = vpt_outputs(vpt_dir, outdir / "staged" / "vpt")
    region_dir = staged_region_dir(
        raw_region_dir, vpt_dir, outdir / "staged" / raw_region_dir.name
    )

    # Elements are named <slide>_<region>_<element>. Both halves are passed rather than
    # defaulted: merscope() takes the slide from the parent directory's name and the region
    # from the directory it read, so defaulting would let the names follow wherever the
    # files sit. Passing them keeps a staged sample's elements identical to a direct read.
    with timer("Read MERSCOPE"):
        sdata = merscope(
            path=region_dir,
            vpt_outputs=paths,
            slide_name=args.sample,
            region_name=raw_region_dir.name,
            z_layers=Z_LAYER,
        )

    print("\nElements:")
    for name in sdata.images:
        print(f"  image   {name}")
    for name in sdata.labels:
        print(f"  labels  {name}")
    for name in sdata.shapes:
        print(f"  shapes  {name}")
    for name in sdata.points:
        print(f"  points  {name}")

    # Record the sample id in the object: every later step and the report read it from
    # there rather than parsing a staged filename. The segmentation is recorded beside it
    # because two stores of the same region are otherwise indistinguishable.
    for name, table in sdata.tables.items():
        table.obs["sample"] = args.sample
        table.uns["z_layer"] = Z_LAYER
        table.uns["segmentation"] = "cellpose"
        table.uns["vpt_path"] = str(vpt_dir.resolve())
        print(f"  table   {name}: {table.n_obs:,} cells x {table.n_vars:,} genes")

    with timer("Write Zarr"):
        sdata.write(output_path, overwrite=True)

    timing_summary(outdir / f"{args.sample}.create_spatialdata_cellpose.timing.tsv")


if __name__ == "__main__":
    main()
