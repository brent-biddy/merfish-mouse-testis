#!/usr/bin/env python3
"""
create_spatialdata.py - Convert a Vizgen MERSCOPE output directory to a SpatialData Zarr store.

Reads a MERSCOPE region directory (cell_by_gene.csv, cell_metadata.csv, the
cell boundaries, the detected transcripts, and the mosaic images) and writes it
out as a SpatialData Zarr store. No cells or genes are filtered here; filtering
happens in later steps, which read the store rather than the MERSCOPE directory.

Boundaries written as a pre-VPT cell_boundaries/ directory of per-FOV HDF5 files
are converted to the cell_boundaries.parquet the reader expects; a region with
neither is an error rather than a store with no polygons in it.

Writes <outdir>/<sample>.create_spatialdata.zarr plus a timing TSV.

Usage:
    create_spatialdata.py --sample testis_01 \\
        --path data/raw/testis_01 \\
        --outdir results/testis_01/create_spatialdata
"""

import argparse
from pathlib import Path

import geopandas as gpd
import h5py
import pandas as pd
from shapely.geometry import MultiPolygon, Polygon
from spatialdata_io import merscope

from timer import timer, timing_summary

# merscope()'s z_layers default, pinned so a library change cannot move the plane silently.
# It selects the image alone: one per-stain mosaic TIFF set, read into the element
# <slide>_<region>_z<n>. Recorded into the table below as provenance; nothing reads it back.
Z_LAYER = 3

# The plane the HDF5 conversion reads, and the ZIndex it stamps on the parquet. Not a
# merscope() argument: its reader keeps only ZIndex 0 -- undocumented, in _get_polygons --
# so any other value converts cleanly, loses every polygon, and dies in pandas.
BOUNDARY_Z_INDEX = 0


def convert_hdf5_boundaries(boundaries_dir, output_path):
    """Write per-FOV HDF5 cell boundaries as the single parquet merscope() reads.

    Pre-VPT MERSCOPE runs write cell_boundaries/feature_data_<fov>.hdf5 instead.
    """
    entity_ids = []
    geometries = []

    for path in sorted(boundaries_dir.glob("feature_data_*.hdf5")):
        with h5py.File(path, "r") as handle:
            for entity_id, cell in handle["featuredata"].items():
                plane = cell.get(f"zIndex_{BOUNDARY_Z_INDEX}")
                if plane is None:
                    continue

                # A cell segmented into several pieces on this plane has p_0, p_1, ..., so
                # parts become one MultiPolygon. Under three vertices cannot close into a
                # ring, and shapely raises on one or two of them; a cell left with no parts
                # at all is skipped.
                parts = []
                for part in plane.values():
                    coordinates = part["coordinates"][0]
                    if len(coordinates) >= 3:
                        parts.append(Polygon(coordinates))
                if not parts:
                    continue

                entity_ids.append(entity_id)
                # MultiPolygon even for a single part: merscope() calls MultiPolygon(x.geoms)
                # on every row, and a bare Polygon has no .geoms.
                geometries.append(MultiPolygon(parts))

    # The HDF5 group names are already strings, so they are left alone rather than cast to
    # the int64 a VPT parquet holds: merscope() indexes on str(EntityID), so both land the
    # same.
    boundaries = gpd.GeoDataFrame(
        {"EntityID": entity_ids, "ZIndex": BOUNDARY_Z_INDEX, "Geometry": geometries},
        geometry="Geometry",
    )
    boundaries.to_parquet(output_path)
    print(f"Converted {len(boundaries):,} cell boundaries at z={BOUNDARY_Z_INDEX}.")


def staged_region_dir(region_dir, staging_dir):
    """Return a region directory merscope() can read, staging a fixed-up one if needed.

    Two things it fixes. The CSVs: merscope() hands both straight to AnnData, which
    requires their indexes to match exactly, and some MERSCOPE versions write them in
    different orders — the load then fails with "Index of obs must match index of X".
    The boundaries: older runs write a cell_boundaries/ directory of per-FOV HDF5 files
    instead of the parquet merscope() reads. Returns the input unchanged when neither
    applies.
    """
    # Only the ids are needed; usecols keeps the whole counts matrix out of memory.
    # dtype={0: str} rather than dtype=str: a scalar dtype does not reach index_col.
    counts_index = pd.read_csv(
        region_dir / "cell_by_gene.csv", index_col=0, usecols=[0], dtype={0: str}
    ).index
    metadata = pd.read_csv(region_dir / "cell_metadata.csv", index_col=0, dtype={0: str})
    reorder = not counts_index.equals(metadata.index)

    parquet = region_dir / "cell_boundaries.parquet"
    hdf5_dir = region_dir / "cell_boundaries"
    if not parquet.exists() and not hdf5_dir.is_dir():
        raise FileNotFoundError(
            f"No cell boundaries in {region_dir}: expected either {parquet.name} or a "
            f"{hdf5_dir.name}/ directory of per-FOV HDF5 files. Without them merscope() "
            f"loads no polygons and writes a store no spatial figure can be drawn from."
        )
    convert = not parquet.exists()

    if not reorder and not convert:
        return region_dir

    # Symlink the rest of the region — the mosaic images are far too large to copy — and
    # write only what is being fixed, so the instrument output is never modified.
    staging_dir.mkdir(parents=True, exist_ok=True)
    replaced = set()
    if reorder:
        replaced.add("cell_metadata.csv")
    if convert:
        replaced.add(hdf5_dir.name)
    for entry in region_dir.iterdir():
        if entry.name not in replaced:
            # Replaced rather than skipped: exists() follows the link, so one left by an
            # earlier run pointing somewhere gone reads as absent and symlink_to then fails.
            link = staging_dir / entry.name
            link.unlink(missing_ok=True)
            link.symlink_to(entry.resolve())

    if reorder:
        print(f"Reordering cell_metadata.csv to match cell_by_gene.csv ({len(counts_index):,} cells).")
        metadata.loc[counts_index].to_csv(staging_dir / "cell_metadata.csv")

    if convert:
        print(f"Converting {hdf5_dir.name}/ to cell_boundaries.parquet.")
        with timer("Convert boundaries"):
            convert_hdf5_boundaries(hdf5_dir, staging_dir / "cell_boundaries.parquet")

    return staging_dir


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a Vizgen MERSCOPE output directory to a SpatialData Zarr store"
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="Sample identifier",
    )
    parser.add_argument(
        "--path",
        required=True,
        help="MERSCOPE region output directory",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory to write <sample>.create_spatialdata.zarr into (default: current directory)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_path = outdir / f"{args.sample}.create_spatialdata.zarr"

    print(f"Sample:  {args.sample}")
    print(f"Input:   {args.path}")
    print(f"Output:  {output_path}")
    print(f"Z-layer: {Z_LAYER}")

    raw_region_dir = Path(args.path)
    region_dir = staged_region_dir(raw_region_dir, outdir / "staged" / raw_region_dir.name)

    # Elements are named <slide>_<region>_<element>. Both halves are passed rather than
    # defaulted: merscope() takes the slide from the parent directory's name and the region
    # from the directory it read, so defaulting would let the names follow wherever the
    # files sit. Passing them keeps a staged sample's elements identical to a direct read.
    with timer("Read MERSCOPE"):
        sdata = merscope(
            path=region_dir,
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
    # there rather than parsing a staged filename.
    for name, table in sdata.tables.items():
        table.obs["sample"] = args.sample
        table.uns["z_layer"] = Z_LAYER
        print(f"  table   {name}: {table.n_obs:,} cells x {table.n_vars:,} genes")

    with timer("Write Zarr"):
        sdata.write(output_path, overwrite=True)

    timing_summary(outdir / f"{args.sample}.create_spatialdata.timing.tsv")


if __name__ == "__main__":
    main()
