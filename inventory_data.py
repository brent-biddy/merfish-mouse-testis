#!/usr/bin/env python3
"""Find the MERSCOPE data under a tree, and list the least of it the pipeline needs.

Stdlib only, so it runs on a login node without the container. Each input is either a directory,
walked for anything the pipeline can read, or a samplesheet, whose path columns are taken as given:

    python3 inventory_data.py /ourdisk/hpc/lilab/babiddy/dont_archive/alex_vizgen
    python3 inventory_data.py assets/samplesheet_cellpose.csv

A directory is recognised by what is in it, not by its name:

    region    images/micron_to_mosaic_pixel_transform.csv    MERSCOPE output, steps 1, 1a, 1b
    vpt       *_micron_space.parquet                          VPT or Vizgen cellpose, step 1b
    cellpose  labels.npy                                      the lab's merged cellpose, step 1a

--steps narrows the files counted as needed to the steps you mean to run; which ones a file serves
is printed beside it. --globus writes those files as a batch for `globus transfer --batch`, each
path relative to --root on both sides, so the tree keeps its shape on the far end.
"""

import argparse
import csv
import os
import shlex
import sys
from pathlib import Path

# The plane merscope() reads images from -- Z_LAYER in create_spatialdata*.py.
Z_LAYER = 3

# kind -> [(label, patterns, steps)]. A label is satisfied when any pattern matches; patterns are
# globs relative to the directory. Step 1b takes only the images and transform from a region: its
# counts, metadata and boundaries come from vpt_outputs, and its transcripts are staged over the
# region's own.
RULES = {
    "region": [
        ("counts", ["cell_by_gene.csv"], {"1", "1a"}),
        ("metadata", ["cell_metadata.csv"], {"1"}),
        ("boundaries", ["cell_boundaries.parquet", "cell_boundaries/feature_data_*.hdf5"], {"1"}),
        ("transcripts", ["detected_transcripts.csv"], {"1", "1a"}),
        ("transform", ["images/micron_to_mosaic_pixel_transform.csv"], {"1", "1a", "1b"}),
        ("mosaic images", [f"images/mosaic_*_z{Z_LAYER}.tif"], {"1", "1a", "1b"}),
    ],
    "vpt": [
        ("counts", ["cellpose_cell_by_gene.csv"], {"1b"}),
        ("metadata", ["cellpose_cell_metadata.csv"], {"1b"}),
        ("boundaries", ["cellpose_micron_space.parquet"], {"1b"}),
        ("transcripts", ["detected_transcripts.csv"], {"1b"}),
    ],
    "cellpose": [
        ("labels", ["labels.npy"], {"1a"}),
        ("label map", ["*-cellp-label-map.npz"], {"1a"}),
        ("counts", ["*-cell-by-gene.npz"], {"1a"}),
        ("bounding boxes", ["*-slicee.pkl"], {"1a"}),
        ("centres", ["cell-coms.npz"], {"1a"}),
        ("volumes", ["cell-vols.npz"], {"1a"}),
    ],
}
# A step 1a output is itself a vpt directory, for step 1b; images reach 1b through the region.
COLUMN_KIND = {"path": "region", "vpt_path": "vpt", "cellpose_path": "cellpose"}
ALL_STEPS = {"1", "1a", "1b"}


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n} B" if unit == "B" else f"{n:,.1f} {unit}"
        n /= 1024


def size(p):
    try:
        return p.stat().st_size
    except OSError:
        return 0


def classify(directory, filenames):
    names = set(filenames)
    if (directory / "images" / "micron_to_mosaic_pixel_transform.csv").is_file():
        return "region"
    if "labels.npy" in names:
        return "cellpose"
    if any(n.endswith("_micron_space.parquet") for n in names):
        return "vpt"
    return None


def discover(root):
    """Yield (kind, directory) for everything under root the pipeline can read."""
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True, onerror=lambda e: print(f"  !! {e}")):
        dirnames.sort()
        kind = classify(Path(dirpath), filenames)
        if kind:
            yield kind, Path(dirpath)
        # Per-FOV and image directories hold thousands of files and never another dataset.
        dirnames[:] = [d for d in dirnames if d not in ("images", "cell_boundaries", "result_tiles")]


def from_samplesheet(sheet):
    with open(sheet, newline="") as fh:
        for row in csv.DictReader(fh):
            for column, kind in COLUMN_KIND.items():
                if row.get(column):
                    yield kind, Path(row[column]), row["sample"]


def inventory(kind, directory, steps):
    """Print one directory's report; return (needed files, bytes read by any step, missing count)."""
    if not directory.is_dir():
        print("    !! not a directory (missing, or a broken link)")
        return [], 0, 1

    needed, matched, missing, read_bytes = [], set(), 0, 0
    for label, patterns, used_by in RULES[kind]:
        hits = sorted({p for pat in patterns for p in directory.glob(pat) if p.is_file()})
        tag = "step " + ", ".join(sorted(used_by))
        wanted = bool(used_by & steps)
        if not hits:
            if wanted:
                missing += 1
            print(f"    {'MISSING' if wanted else 'absent':<8} {label:<15} {' or '.join(patterns):<50} {'':>10}  {tag}")
            continue
        total = sum(size(p) for p in hits)
        read_bytes += total
        matched.update(hits)
        if wanted:
            needed.extend(hits)
        shown = hits[0].relative_to(directory) if len(hits) == 1 else f"{len(hits)} files"
        print(f"    {'need' if wanted else 'skip':<8} {label:<15} {str(shown):<50} {human(total):>10}  {tag}")

    unread = {}
    for p in directory.rglob("*"):
        if p.is_file() and p not in matched:
            top = p.relative_to(directory).parts[0]
            unread[top] = unread.get(top, 0) + size(p)
    if unread:
        ranked = sorted(unread.items(), key=lambda kv: -kv[1])
        print(f"    not read by any step, {human(sum(unread.values()))}:")
        for top, n in ranked[:8]:
            print(f"      {human(n):>10}  {top}")
        if len(ranked) > 8:
            print(f"      {human(sum(n for _, n in ranked[8:])):>10}  ... {len(ranked) - 8} more entries")
    print(f"    needed: {human(sum(size(p) for p in needed))} of {human(read_bytes)} the pipeline can read")
    return needed, read_bytes, missing


def write_globus(path, files, root):
    with open(path, "w") as fh:
        for f in files:
            rel = shlex.quote(str(f.relative_to(root)))
            fh.write(f"{rel} {rel}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs", nargs="+", type=Path, help="directories to walk, or samplesheets")
    parser.add_argument("--steps", default="1,1a,1b",
                        help="comma-separated steps to keep files for (default: %(default)s)")
    parser.add_argument("--globus", type=Path, help="write the needed files here as a globus --batch file")
    parser.add_argument("--root", type=Path,
                        help="base path the batch lines are relative to (default: the inputs' common path)")
    args = parser.parse_args()

    steps = set(args.steps.split(","))
    if not steps <= ALL_STEPS:
        parser.error(f"--steps takes {', '.join(sorted(ALL_STEPS))}, got {args.steps}")

    found = []
    for item in args.inputs:
        if item.suffix == ".csv":
            found.extend(from_samplesheet(item))
        else:
            found.extend((kind, d, None) for kind, d in discover(item))

    needed, seen, readable, missing, by_kind = [], set(), 0, 0, {}
    for kind, directory, sample in found:
        key = (kind, directory.resolve())
        if key in seen:
            continue
        seen.add(key)
        print(f"\n[{kind}] {directory}" + (f"  (sample {sample})" if sample else ""))
        files, n, m = inventory(kind, directory, steps)
        needed.extend(files)
        readable += n
        missing += m
        by_kind.setdefault(kind, [0, 0])
        by_kind[kind][0] += 1
        by_kind[kind][1] += sum(size(p) for p in files)

    # A file reached by two routes -- a region named in two samplesheets -- is sent once.
    needed = sorted({p.absolute() for p in needed})
    total = sum(size(p) for p in needed)

    print(f"\n== Summary for steps {', '.join(sorted(steps))}")
    for kind, (count, n) in sorted(by_kind.items()):
        print(f"  {count:>3} {kind:<9} directories  {human(n):>10} needed")
    print(f"  {len(needed):,} files, {human(total)} to transfer ({human(readable)} readable by any step)")
    if missing:
        print(f"  {missing} needed file(s) missing -- see MISSING above.")

    links = [p for p in needed if p.is_symlink()]
    if links:
        print(f"  {len(links)} needed file(s) are symlinks; globus does not follow them. First: {links[0]}")

    if args.globus:
        root = (args.root or Path(os.path.commonpath([str(p) for p in needed]))).absolute()
        if root.is_file():
            root = root.parent
        write_globus(args.globus, needed, root)
        print(f"\n  wrote {args.globus}, paths relative to {root}")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
