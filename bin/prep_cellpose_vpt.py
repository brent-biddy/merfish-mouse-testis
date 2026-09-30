#!/usr/bin/env python3
"""
prep_cellpose_vpt.py - Convert a bespoke cellpose segmentation of a MERSCOPE region to
the three VPT output files create_spatialdata_cellpose.py reads.

Not a VPT run: this is the lab's own cellpose pipeline (area_cellpose.py per FOV, merge.py
to stitch), which writes numpy arrays where merscope() wants CSVs and a parquet. Writing
those here keeps one reader downstream.

Inputs, all in the merged segmentation directory:

    labels.npy                  (7, height, width) uint32 raw labels, full mosaic pixels
      or <sample>-full-merged-labels-compressed.npz, the same array zipped, ~1% the size
    <sample>-cellp-label-map.npz    raw label -> compact cell id, indexed by label - 1
    <sample>-cell-by-gene.npz   counts, row 0 is background
    <sample>-slicee.pkl         per-label bounding boxes, merge.py's find_objects output
    cell-coms.npz               (z, y, x) centre per cell, mosaic pixels
    cell-vols.npz               voxels per cell

None of the index conventions are recoverable from the files, so, from merge.py:

    compact_id = label_map[raw_label - 1]    0 means the label is absent
    counts     = cbg[compact_id]             row 0 is background, not a cell
    centre     = cell_coms[compact_id - 1]
    volume     = cell_vols[compact_id - 1]

label_map is indexed by `label - 1` because merge.py builds it over ndi.find_objects
output, whose element i describes label i + 1. Indexing by the label instead is off by one
wherever a label is absent, silently: it gives a cell its neighbour's counts. check_volumes
is what catches that before anything is written.

The segmentation is 3D over 7 planes and a store holds one, so each cell's boundary comes
from the plane it covers most; a cell centred near the top or bottom of the section still
gets one, and which plane lands in obs["z_plane"].

Everything that needs the raster is done in one pass over its planes -- the areas, the
transcript lookup and each cell's mask -- because the compressed copy can only be read in
order: a second pass would mean inflating ~45 GB again. labels.npy is used when present, as
it is faster to read; otherwise the compressed copy, so a transfer can leave labels.npy out.

Every cell the segmentation found gets a row in the counts and the metadata, including one
whose mask traced to nothing: the parquet is the only output a missing polygon can keep it
out of, and metadata["has_boundary"] says which those are.

Writes <outdir>/cellpose_{cell_by_gene.csv,cell_metadata.csv,micron_space.parquet}, a
detected_transcripts.csv carrying the cell_id column VPT would have written, and a timing
TSV. Point create_spatialdata_cellpose.py --vpt_path at <outdir>.

Usage:
    prep_cellpose_vpt.py --sample b2r0_cellpose \\
        --path data/raw/MsTestis_B2/region_0 \\
        --cellpose_path data/cellpose/b2r0-051223 \\
        --outdir results/b2r0_cellpose/prep_cellpose_vpt
"""

import argparse
import pickle
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import MultiPolygon, Polygon
from skimage.measure import find_contours

from timer import timer, timing_summary

# merscope() keeps only ZIndex 0 out of whatever the boundary file holds, so every polygon
# is written at 0 whichever plane it was taken from.
BOUNDARY_Z_INDEX = 0

# Mosaic z spacing. The MERSCOPE records the pixel size but not this, so it cannot be read
# from the region; merge.py's anisotropy constant assumes the same. Only the volume column
# depends on it: wrong here, volumes scale and nothing else moves.
Z_STEP_MICRONS = 1.5

# The names create_spatialdata_cellpose.py looks for: a Vizgen cellpose delivery's, where
# all three carry the method. A stock VPT run prefixes only the boundaries.
OUTPUT_NAMES = {
    "cell_by_gene": "cellpose_cell_by_gene.csv",
    "cell_metadata": "cellpose_cell_metadata.csv",
    "cell_boundaries": "cellpose_micron_space.parquet",
}

# Written by merge.py under the sample's own prefix, which is not the sample id this
# pipeline uses -- matched by suffix instead of named.
LABEL_MAP_SUFFIX = "-cellp-label-map.npz"
COUNTS_SUFFIX = "-cell-by-gene.npz"
COMPRESSED_LABELS_SUFFIX = "-full-merged-labels-compressed.npz"

# Read from the region and rewritten with a cell_id column. VPT's own name for both the
# column and the file, and -1 for a transcript in no cell -- see vpt/partition_transcripts.
TRANSCRIPTS_NAME = "detected_transcripts.csv"
UNASSIGNED = -1


def find_one(directory, suffix):
    """Return the single file in directory ending in suffix."""
    matches = sorted(directory.glob(f"*{suffix}"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly one *{suffix} in {directory}, found {len(matches)}. "
            f"merge.py writes one per merged segmentation."
        )
    return matches[0]


class LabelPlanes:
    """The label raster, one z plane at a time, from labels.npy or merge.py's compressed copy.

    labels.npy is memory-mapped. The compressed copy is a zip holding that same .npy, which
    cannot be mapped -- np.load would inflate every plane at once, ~45 GB for a region -- so
    it is streamed: the array is C-ordered with z first, so each plane is the next
    height * width values. Iterating yields (z, plane); a streamed plane reuses one buffer,
    so it is valid only until the next.
    """

    def __init__(self, cellpose_dir):
        self.path = cellpose_dir / "labels.npy"
        if self.path.exists():
            self.array = np.load(self.path, mmap_mode="r")
            self.shape, self.dtype = self.array.shape, self.array.dtype
            return

        self.path = find_one(cellpose_dir, COMPRESSED_LABELS_SUFFIX)
        self.array = None
        with zipfile.ZipFile(self.path) as archive:
            members = archive.namelist()
            if len(members) != 1:
                raise ValueError(f"{self.path.name} holds {members}; expected one array.")
            self.member = members[0]
            with archive.open(self.member) as stream:
                self.shape, fortran_order, self.dtype = self._read_header(stream)
        if fortran_order or len(self.shape) != 3:
            raise ValueError(
                f"{self.path.name} holds a {'Fortran-ordered ' if fortran_order else ''}"
                f"{self.shape} array; streaming needs a C-ordered (z, height, width) one."
            )

    @staticmethod
    def _read_header(stream):
        version = np.lib.format.read_magic(stream)
        if version == (1, 0):
            return np.lib.format.read_array_header_1_0(stream)
        return np.lib.format.read_array_header_2_0(stream)

    def __iter__(self):
        if self.array is not None:
            for z in range(self.shape[0]):
                yield z, np.asarray(self.array[z])
            return

        plane = np.empty(self.shape[1:], dtype=self.dtype)
        view = memoryview(plane).cast("B")
        with zipfile.ZipFile(self.path) as archive, archive.open(self.member) as stream:
            self._read_header(stream)
            for z in range(self.shape[0]):
                filled = 0
                while filled < len(view):
                    n = stream.readinto(view[filled:])
                    if not n:
                        raise EOFError(f"{self.path.name} ends partway through plane {z}.")
                    filled += n
                yield z, plane
            # Reading to the end is what makes zipfile check the CRC, so a corrupt copy
            # fails here rather than passing through as wrong labels.
            if stream.read():
                raise ValueError(f"{self.path.name} holds more than {self.shape}.")


def load_segmentation(cellpose_dir):
    """Read the merged segmentation arrays, checking they agree on how many cells there are."""
    labels = LabelPlanes(cellpose_dir)
    label_map = np.load(find_one(cellpose_dir, LABEL_MAP_SUFFIX))["label_map"]
    counts = np.load(find_one(cellpose_dir, COUNTS_SUFFIX))["cbg"]
    centres = np.load(cellpose_dir / "cell-coms.npz")["cell_coms"]
    volumes = np.load(cellpose_dir / "cell-vols.npz")["cell_vols"]

    n_cells = int(label_map.max())
    print(f"Labels:  {labels.shape} {labels.dtype}, from {labels.path.name}")
    print(f"Cells:   {n_cells:,}")

    if not (len(centres) == len(volumes) == n_cells):
        raise ValueError(
            f"Cell counts disagree: label_map.max()={n_cells:,}, "
            f"cell-coms={len(centres):,}, cell-vols={len(volumes):,}."
        )
    if len(counts) != n_cells + 1:
        raise ValueError(
            f"{find_one(cellpose_dir, COUNTS_SUFFIX).name} has {len(counts):,} rows; "
            f"expected {n_cells + 1:,} -- one per cell plus row 0 for background."
        )

    # Every array below is indexed by the compact id, so two labels sharing one would
    # silently merge two cells into a single row.
    present = np.nonzero(label_map)[0]
    if len(np.unique(label_map[present])) != len(present):
        raise ValueError(
            f"label_map is not one-to-one: {len(present):,} labels map to "
            f"{len(np.unique(label_map[present])):,} distinct cells."
        )
    return labels, label_map, counts, centres, volumes


def label_of_cell(label_map):
    """Return an array giving each compact cell id the label_map index that reaches it.

    The inverse of the label_map, which load_segmentation has checked is one-to-one. Index
    0 is unused, so the array lines up with the compact ids.
    """
    present = np.nonzero(label_map)[0]
    reverse = np.zeros(int(label_map.max()) + 1, dtype=np.int64)
    reverse[label_map[present]] = present
    return reverse


def cell_bounds(slices, label_index, transform):
    """Return per-cell (min_x, min_y, max_x, max_y) in microns.

    From merge.py's bounding boxes rather than the traced polygons: they cover every plane
    the cell occupies, and they exist for a cell whose trace produced nothing.
    """
    row_starts = np.array([slices[index][1].start for index in label_index])
    row_stops = np.array([slices[index][1].stop for index in label_index])
    column_starts = np.array([slices[index][2].start for index in label_index])
    column_stops = np.array([slices[index][2].stop for index in label_index])

    # Both corners through the affine, then min and max, so a rotation could not flip them.
    x0, y0 = to_microns(row_starts, column_starts, transform)
    x1, y1 = to_microns(row_stops, column_stops, transform)
    return (np.minimum(x0, x1), np.minimum(y0, y1),
            np.maximum(x0, x1), np.maximum(y0, y1))


def transcript_pixels(transcripts, to_pixels, shape):
    """Return each transcript's (rows, columns, planes) in the raster, and which are on it.

    global_z is a plane index, 0 to 6, not a micron depth -- so it selects the labels plane
    directly, and the transcript's micron position indexes into it.
    """
    pixels = to_pixels @ np.stack([
        transcripts["global_x"].to_numpy(),
        transcripts["global_y"].to_numpy(),
        np.ones(len(transcripts)),
    ])
    columns = np.rint(pixels[0]).astype(np.int64)
    rows = np.rint(pixels[1]).astype(np.int64)
    planes = transcripts["global_z"].to_numpy().astype(np.int64)

    n_planes, height, width = shape
    on_raster = (
        (rows >= 0) & (rows < height) & (columns >= 0) & (columns < width)
        & (planes >= 0) & (planes < n_planes)
    )
    return rows, columns, planes, on_raster


def scan_planes(labels, label_map, slices, pixels):
    """Do everything that needs the raster in one pass over its planes.

    Returns
      areas    (n planes, n labels + 1) pixels per raw label per plane
      cell_id  per transcript, the cell owning the voxel it lands in, or UNASSIGNED
      masks    per label - 1, the cell's mask on the plane it covers most, bit-packed

    A transcript is in whatever cell owns its voxel, the same raster merge.py counted, so
    the totals should reproduce its cell-by-gene matrix. Assignment is 3D and a boundary is
    one plane, so a transcript can belong to a cell and still fall outside its polygon.

    Which plane a cell covers most is only known after the last one, so each cell keeps its
    largest mask so far and a later plane replaces it only by covering strictly more -- a
    tie keeps the earlier plane, the one argmax over the areas picks.
    """
    n_labels = len(label_map)
    areas = np.zeros((labels.shape[0], n_labels + 1), dtype=np.int64)
    largest = np.zeros(n_labels, dtype=np.int64)     # per label - 1, area of the kept mask
    masks = {}

    rows, columns, planes, on_raster = pixels
    cell_id = np.full(len(rows), UNASSIGNED, dtype=np.int64)

    for z, plane in labels:
        # bincount takes intp, so it casts the whole plane -- twice the plane again on top
        # of the plane itself. In blocks that cast is a sixteenth, for identical counts.
        for block in np.array_split(plane.ravel(), 16):
            areas[z] += np.bincount(block, minlength=n_labels + 1)

        here = np.nonzero(on_raster & (planes == z))[0]
        raw = plane[rows[here], columns[here]]
        in_cell = raw > 0
        # label_map holds 0 for a label no cell claims, which is not VPT's sentinel.
        mapped = label_map[raw[in_cell] - 1]
        found = here[in_cell][mapped > 0]
        cell_id[found] = mapped[mapped > 0]

        grown = np.nonzero((areas[z][1:] > largest) & (label_map > 0))[0]
        for index in grown:
            _, cell_rows, cell_columns = slices[index]
            mask = plane[cell_rows, cell_columns] == index + 1
            masks[index] = (np.packbits(mask), mask.shape)
        largest[grown] = areas[z][grown + 1]

        print(f"  z={z}: {int((areas[z][1:] > 0).sum()):,} labels present, "
              f"{len(found):,} of {len(here):,} transcripts in a cell, "
              f"{len(grown):,} masks kept")
    return areas, cell_id, masks


def check_volumes(areas, label_map, volumes):
    """Fail unless each label's voxels match the volume recorded for the cell it maps to."""
    present = np.nonzero(label_map)[0]          # index into label_map, i.e. label - 1
    from_raster = areas[:, present + 1].sum(axis=0)
    recorded = volumes[label_map[present] - 1]

    mismatches = int((from_raster != recorded).sum())
    if mismatches:
        raise ValueError(
            f"{mismatches:,} of {len(present):,} cells have a raster volume that differs "
            f"from cell-vols.npz. The label map, the raster and the volumes are not from "
            f"one merge.py run, or its indexing has changed."
        )
    print(f"Volume check: {len(present):,} cells match cell-vols.npz exactly.")


def micron_to_pixel_transform(region_dir):
    """Return the region's own affine taking microns to mosaic pixels."""
    path = region_dir / "images" / "micron_to_mosaic_pixel_transform.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"No {path.name} in {path.parent}. It is what relates the segmentation's pixel "
            f"coordinates to the micron space the transcripts and images use."
        )
    return pd.read_csv(path, sep=r"\s+", header=None).values


def check_counts(transcripts, counts, genes, entity_ids):
    """Report how far the transcripts' own per-cell totals sit from merge.py's matrix.

    Printed rather than raised: merge.py's assignment rule is not recorded anywhere, and
    this has never run against a real segmentation. Tighten it to an error once a run has
    shown the two agree exactly.
    """
    assigned = transcripts[transcripts["cell_id"] > 0]
    observed = (
        pd.crosstab(assigned["cell_id"], assigned["gene"])
        .reindex(index=entity_ids, columns=genes, fill_value=0)
        .to_numpy()
    )
    recorded = counts[entity_ids]

    per_cell = (observed == recorded).all(axis=1)
    difference = int(np.abs(observed - recorded).sum())
    print(f"Counts check: {int(per_cell.sum()):,} of {len(entity_ids):,} cells reproduce "
          f"cell-by-gene exactly; {difference:,} transcripts differ in total.")
    if difference:
        print("  merge.py assigned transcripts by some other rule -- the counts written "
              "here are still its own, so only the cell_id column is in question.")


def to_microns(rows, columns, transform):
    """Map mosaic pixel (row, column) arrays to micron (x, y) arrays."""
    pixels = np.stack([columns, rows, np.ones_like(columns)])
    x, y, _ = transform @ pixels
    return x, y


def cell_polygon(mask, y_offset, x_offset, transform):
    """Return the MultiPolygon of a boolean cell mask, in microns, or None.

    Padded before tracing so a cell touching its bounding box still closes. Several pieces
    on one plane trace as several contours; under three vertices is not a ring.
    """
    contours = find_contours(np.pad(mask, 1), 0.5)

    parts = []
    for contour in contours:
        if len(contour) < 3:
            continue
        rows = contour[:, 0] - 1 + y_offset
        columns = contour[:, 1] - 1 + x_offset
        x, y = to_microns(rows, columns, transform)
        parts.append(Polygon(np.column_stack([x, y])))

    if not parts:
        return None

    # A hole traces as its own contour, so a cell with one comes back as nested shells:
    # invalid, and about half again too large. make_valid turns those into interior rings.
    shape = shapely.make_valid(MultiPolygon(parts))

    # MultiPolygon even for one piece: merscope() calls MultiPolygon(x.geoms) on every row,
    # and a bare Polygon has no .geoms. make_valid returns one whenever the cell is single.
    polygons = [part for part in getattr(shape, "geoms", [shape])
                if part.geom_type == "Polygon"]
    if not polygons:
        return None
    return MultiPolygon(polygons)


def trace_boundaries(masks, slices, label_map, transform):
    """Trace one polygon per cell from the mask scan_planes kept for it.

    A cell whose mask traces to nothing is left out of the parquet; it still gets a row in
    the counts and the metadata, so the table carries every cell the segmentation found.
    """
    entity_ids = []
    geometries = []

    for index in sorted(masks):
        packed, shape = masks[index]
        mask = np.unpackbits(packed, count=shape[0] * shape[1]).reshape(shape).astype(bool)
        _, rows, columns = slices[index]
        polygon = cell_polygon(mask, rows.start, columns.start, transform)
        if polygon is None:
            continue

        entity_ids.append(int(label_map[index]))
        geometries.append(polygon)

    return gpd.GeoDataFrame(
        {"EntityID": entity_ids, "ZIndex": BOUNDARY_Z_INDEX, "Geometry": geometries},
        geometry="Geometry",
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a bespoke cellpose segmentation of a MERSCOPE region to VPT "
                    "output files"
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="Sample identifier",
    )
    parser.add_argument(
        "--path",
        required=True,
        help="MERSCOPE region output directory, for the gene names and the micron transform",
    )
    parser.add_argument(
        "--cellpose_path",
        required=True,
        help="Merged cellpose segmentation directory, holding labels.npy (or its compressed "
             "copy) and merge.py's arrays",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory to write the VPT files into (default: current directory)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    region_dir = Path(args.path)
    cellpose_dir = Path(args.cellpose_path)

    print(f"Sample:   {args.sample}")
    print(f"Region:   {region_dir}")
    print(f"Cellpose: {cellpose_dir}")
    print(f"Output:   {outdir}")

    with timer("Load segmentation"):
        labels, label_map, counts, centres, volumes = load_segmentation(cellpose_dir)

    # Gene order is the region's own, since the counts array carries no column names.
    genes = pd.read_csv(region_dir / "cell_by_gene.csv", index_col=0, nrows=0).columns
    if len(genes) != counts.shape[1]:
        raise ValueError(
            f"{region_dir / 'cell_by_gene.csv'} names {len(genes):,} genes but the "
            f"segmentation counts have {counts.shape[1]:,} columns."
        )
    print(f"Genes:    {len(genes):,}")

    # Both directions are wanted: transcripts come in microns and index the pixel raster,
    # while the segmentation is in pixels and everything written out is in microns.
    to_pixels = micron_to_pixel_transform(region_dir)
    transform = np.linalg.inv(to_pixels)

    # merge.py already ran find_objects over the volume and pickled it; reuse, do not rebuild.
    slicee_path = find_one(cellpose_dir, "-slicee.pkl")
    with open(slicee_path, "rb") as handle:
        slices = pickle.load(handle)
    if len(slices) != len(label_map):
        raise ValueError(
            f"{slicee_path.name} has {len(slices):,} entries and the label map "
            f"{len(label_map):,}; both index raw label - 1 and must agree."
        )

    with timer("Read transcripts"):
        transcripts = pd.read_csv(region_dir / TRANSCRIPTS_NAME)
        pixels = transcript_pixels(transcripts, to_pixels, labels.shape)

    with timer("Scan planes"):
        areas, cell_id, masks = scan_planes(labels, label_map, slices, pixels)
    transcripts["cell_id"] = cell_id

    check_volumes(areas, label_map, volumes)

    # A label only exists because it has voxels, so none is absent everywhere -- but argmax
    # would quietly answer 0 if one were.
    per_label = areas[:, 1:]
    best_plane = per_label.argmax(axis=0)
    best_plane[label_map == 0] = -1

    with timer("Trace boundaries"):
        boundaries = trace_boundaries(masks, slices, label_map, transform)
    del masks

    # merscope() drops invalid geometries without warning, so a polygon that survives the
    # trace but not the reader would leave its cell in the table with nothing to draw.
    invalid = int((~boundaries.geometry.is_valid).sum())
    if invalid:
        raise ValueError(
            f"{invalid:,} of {len(boundaries):,} traced polygons are invalid. merscope() "
            f"would drop them silently and leave their counts in the table."
        )

    # Every cell the segmentation found, whether or not its mask traced: the counts and the
    # centre do not depend on the polygon, and a cell missing from the table cannot be
    # recovered downstream, where one missing from the parquet is only undrawable.
    entity_ids = np.unique(label_map[label_map > 0])
    label_index = label_of_cell(label_map)[entity_ids]
    order = np.argsort(boundaries["EntityID"].to_numpy())
    print(f"Traced {len(boundaries):,} boundaries for {len(entity_ids):,} cells.")

    with timer("Write counts"):
        cell_by_gene = pd.DataFrame(counts[entity_ids], index=entity_ids, columns=genes)
        cell_by_gene.index.name = "cell"
        cell_by_gene.to_csv(outdir / OUTPUT_NAMES["cell_by_gene"])

    with timer("Write metadata"):
        centre_x, centre_y = to_microns(
            centres[entity_ids - 1, 1], centres[entity_ids - 1, 2], transform
        )
        min_x, min_y, max_x, max_y = cell_bounds(slices, label_index, transform)
        # merge.py counted voxels; a Vizgen cell_metadata.csv holds cubic microns, which is
        # what everything downstream reads this column expecting.
        voxel_microns = abs(np.linalg.det(transform[:2, :2])) * Z_STEP_MICRONS
        metadata = pd.DataFrame(
            {
                "EntityID": entity_ids,
                "volume": volumes[entity_ids - 1] * voxel_microns,
                "volume_voxels": volumes[entity_ids - 1],
                "center_x": centre_x,
                "center_y": centre_y,
                "min_x": min_x,
                "min_y": min_y,
                "max_x": max_x,
                "max_y": max_y,
                # The plane the boundary was taken from, or would have been had it traced;
                # the parquet itself is all ZIndex 0.
                "z_plane": best_plane[label_index],
                "has_boundary": np.isin(entity_ids, boundaries["EntityID"].to_numpy()),
            }
        ).set_index("EntityID")
        metadata.to_csv(outdir / OUTPUT_NAMES["cell_metadata"])

    with timer("Write boundaries"):
        boundaries.iloc[order].to_parquet(outdir / OUTPUT_NAMES["cell_boundaries"])

    with timer("Write transcripts"):
        transcripts.to_csv(outdir / TRANSCRIPTS_NAME, index=False)

    check_counts(transcripts, counts, genes, entity_ids)

    print(f"\nWrote {len(entity_ids):,} cells, "
          f"{len(entity_ids) - len(boundaries):,} of them without a boundary:")
    for name in list(OUTPUT_NAMES.values()) + [TRANSCRIPTS_NAME]:
        print(f"  {name}")

    timing_summary(outdir / f"{args.sample}.prep_cellpose_vpt.timing.tsv")


if __name__ == "__main__":
    main()
