#!/usr/bin/env python3
"""
split_cellpose_section.py - Cut one merge.py segmentation that holds two tissue sections
into one segmentation directory per section, in the layout prep_cellpose_vpt.py reads.

a1r0 is the case: two sections sit in one region, so merge.py wrote one raster, one label
map and one counts matrix for both. prep_cellpose_vpt.py takes a region as one sample, so
the segmentation is split here, before anything is processed.

Inputs, all in --cellpose_path, named by --prefix (a1r0):

    labels.npy                  (7, height, width) uint32 raw labels
    <prefix>-cellp-label-map.npz    raw label -> compact cell id, indexed by label - 1
    <prefix>-cell-by-gene.npz   counts, row 0 is background
    <prefix>-slicee.pkl         per-label bounding boxes
    cell-coms.npz               (z, y, x) centre per cell, mosaic pixels
    cell-vols.npz               voxels per cell

The cut is --split_y, a global_y in microns (Alex Vargo's xscript-metafile-creation.py puts
a1r0's at 3520). Transcripts follow his rule exactly: the lower section is global_y >
split_y and the upper is the rest. Cells follow the same line in pixels: the region's
transform takes split_y to a mosaic row (33500 for a1r0), a cell belongs to the upper
section when its centre is above it, and the cut in the raster is that row less the raster's
row offset. No cell may cross it: a bounding box that does stops the run.

Each section keeps the raw labels but renumbers its cells 1..n, so the label map, counts,
centres and volumes are cut to its cells and nothing else changes. The lower section's
raster starts at the cut, so its rows shift and so does its label offset; the upper's is
unchanged. The offsets to put in the samplesheet are printed and written to
split_cellpose_section.json.

Writes, for each of --upper_name and --lower_name:

    <outdir>/<name>/            the files above, with <name> in place of <prefix>
    <outdir>/regions/<name>/    a copy of --region_path whose detected_transcripts.csv holds
                                only that section's transcripts; everything else is linked

so prep_cellpose_vpt.py takes the first as --cellpose_path and the second as --path and
needs no change. Only labels.npy and the transcripts are large: the raster is copied in
bands, so memory stays small but the disk needs the whole raster again, and the transcripts
are read whole, a few GB.

Usage:
    split_cellpose_section.py --cellpose_path cellposetest --prefix a1r0 \\
        --region_path A1/region_0 --label_offset=-245,-290 --split_y 3520 \\
        --upper_name a1r0-upper --lower_name a1r0-lower --outdir a1r0_split
"""

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

TRANSCRIPTS_NAME = "detected_transcripts.csv"

# Rows of the raster copied at a time: 7 planes x 40,000 columns x 4 bytes is ~1 GB.
BAND_ROWS = 1000


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--cellpose_path", required=True, help="merge.py output directory")
    parser.add_argument("--prefix", required=True, help="file prefix, e.g. a1r0")
    parser.add_argument(
        "--label_offset", required=True,
        help="ROWS,COLUMNS of the raster's (0, 0) in the mosaic; write --label_offset=-245,-290",
    )
    parser.add_argument("--region_path", required=True, help="MERSCOPE region directory")
    parser.add_argument("--split_y", type=float, required=True, help="global_y of the cut, microns")
    parser.add_argument("--upper_name", required=True)
    parser.add_argument("--lower_name", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def write_section(name, cells, source, outdir, row_start, row_stop, shift):
    """Write one section's segmentation directory and return its cell count.

    cells are compact ids in the combined segmentation. row_start and row_stop are the
    section's raster rows in the combined raster, and shift is what its rows lose.
    """
    section_dir = outdir / name
    section_dir.mkdir(parents=True, exist_ok=True)

    # Compact ids 1..n in the combined order, and the same raw labels pointing at them.
    n_cells = len(cells)
    label_index = source["label_index"][cells - 1]
    label_map = np.zeros(len(source["label_map"]), dtype=source["label_map"].dtype)
    label_map[label_index] = np.arange(1, n_cells + 1)

    slicee = [None] * len(source["slicee"])
    for index in label_index:
        planes, rows, columns = source["slicee"][index]
        slicee[index] = (planes, slice(rows.start - shift, rows.stop - shift), columns)

    counts = np.vstack([source["counts"][:1], source["counts"][cells]])

    np.savez(section_dir / f"{name}-cellp-label-map.npz", label_map=label_map)
    np.savez(section_dir / f"{name}-cell-by-gene.npz", cbg=counts)
    np.savez(section_dir / "cell-coms.npz", cell_coms=source["coms"][cells - 1])
    np.savez(section_dir / "cell-vols.npz", cell_vols=source["vols"][cells - 1])
    with open(section_dir / f"{name}-slicee.pkl", "wb") as handle:
        pickle.dump(slicee, handle)

    labels = source["labels"]
    n_planes, _, width = labels.shape
    section_labels = np.lib.format.open_memmap(
        section_dir / "labels.npy", mode="w+", dtype=labels.dtype,
        shape=(n_planes, row_stop - row_start, width),
    )
    for band_start in range(row_start, row_stop, BAND_ROWS):
        band_stop = min(band_start + BAND_ROWS, row_stop)
        section_labels[:, band_start - row_start:band_stop - row_start] = labels[:, band_start:band_stop]
    section_labels.flush()

    return n_cells


def write_region(name, region_dir, outdir, transcripts):
    """Write a copy of the region directory with these transcripts, linking the rest.

    The mosaic images are far too large to copy, and nothing else is changed.
    """
    section_region = outdir / "regions" / name
    section_region.mkdir(parents=True, exist_ok=True)

    for entry in region_dir.iterdir():
        if entry.name == TRANSCRIPTS_NAME:
            continue
        link = section_region / entry.name
        if link.is_symlink():
            link.unlink()
        link.symlink_to(entry.resolve())

    transcripts.to_csv(section_region / TRANSCRIPTS_NAME)


def main():
    args = parse_args()
    seg = Path(args.cellpose_path)
    region_dir = Path(args.region_path)
    outdir = Path(args.outdir)
    offset_rows, offset_columns = (int(value) for value in args.label_offset.split(","))

    # Row 1 of the transform takes microns to a mosaic row: row = global_y * slope + intercept.
    transform = np.loadtxt(region_dir / "images" / "micron_to_mosaic_pixel_transform.csv")
    split_row = int(round(args.split_y * transform[1, 1] + transform[1, 2]))
    cut = split_row - offset_rows

    with open(seg / f"{args.prefix}-slicee.pkl", "rb") as handle:
        slicee = pickle.load(handle)
    label_map = np.load(seg / f"{args.prefix}-cellp-label-map.npz")["label_map"]
    counts = np.load(seg / f"{args.prefix}-cell-by-gene.npz")["cbg"]
    coms = np.load(seg / "cell-coms.npz")["cell_coms"]
    vols = np.load(seg / "cell-vols.npz")["cell_vols"]
    labels = np.load(seg / "labels.npy", mmap_mode="r")

    # label_index[c - 1] is the label - 1 that compact id c came from.
    label_index = np.zeros(int(label_map.max()), dtype=np.int64)
    present = np.nonzero(label_map)[0]
    label_index[label_map[present] - 1] = present

    cells = np.arange(1, len(label_index) + 1)
    is_upper = coms[:, 1] < split_row
    upper_cells = cells[is_upper]
    lower_cells = cells[~is_upper]

    upper_stops = max(slicee[index][1].stop for index in label_index[upper_cells - 1])
    lower_starts = min(slicee[index][1].start for index in label_index[lower_cells - 1])
    if upper_stops > cut or lower_starts < cut:
        raise ValueError(
            f"A cell crosses raster row {cut}: upper cells reach row {upper_stops}, lower "
            f"cells start at row {lower_starts}. Check --split_y and --label_offset."
        )

    source = {
        "label_index": label_index, "label_map": label_map, "slicee": slicee,
        "counts": counts, "coms": coms, "vols": vols, "labels": labels,
    }
    n_upper = write_section(
        args.upper_name, upper_cells, source, outdir, 0, cut, 0,
    )
    n_lower = write_section(
        args.lower_name, lower_cells, source, outdir, cut, labels.shape[1], cut,
    )

    # The first column is VPT's unnamed row index; reading and writing it as the index keeps
    # the header exactly as the instrument wrote it.
    transcripts = pd.read_csv(region_dir / TRANSCRIPTS_NAME, index_col=0)
    is_lower = transcripts["global_y"] > args.split_y
    write_region(args.upper_name, region_dir, outdir, transcripts[~is_lower])
    write_region(args.lower_name, region_dir, outdir, transcripts[is_lower])
    n_transcripts = {
        args.upper_name: int((~is_lower).sum()),
        args.lower_name: int(is_lower.sum()),
    }

    offsets = {
        args.upper_name: [offset_rows, offset_columns],
        args.lower_name: [offset_rows + cut, offset_columns],
    }
    record = {
        "cellpose_path": str(seg.resolve()),
        "region_path": str(region_dir.resolve()),
        "split_y": args.split_y,
        "split_row": split_row,
        "raster_cut_row": cut,
        "cells": {args.upper_name: n_upper, args.lower_name: n_lower},
        "transcripts": n_transcripts,
        "label_offset": offsets,
    }
    with open(outdir / "split_cellpose_section.json", "w") as handle:
        json.dump(record, handle, indent=2)

    for name, offset in offsets.items():
        print(f'{name}: {record["cells"][name]:,} cells, {n_transcripts[name]:,} transcripts, '
              f'label_offset "{offset[0]},{offset[1]}"')


if __name__ == "__main__":
    main()
