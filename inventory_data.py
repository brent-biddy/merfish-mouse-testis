#!/usr/bin/env python3
"""Find the MERSCOPE data under a tree, and list the least of it the pipeline needs.

Stdlib only, so it runs on a login node without the container. Each input is a directory, walked
for anything the pipeline can read; a Globus collection path, walked the same way over the globus
CLI; or a samplesheet, whose path columns are taken as given:

    python3 inventory_data.py /ourdisk/hpc/lilab/babiddy/dont_archive/alex_vizgen
    python3 inventory_data.py 13841aad-d64b-4420-a9b3-2e60762b1972:/segmentation
    python3 inventory_data.py assets/samplesheet_cellpose.csv

A collection path needs the globus CLI on PATH and logged in -- `conda activate globus`.

A directory is recognised by what is in it, not by its name:

    region    images/micron_to_mosaic_pixel_transform.csv    MERSCOPE output, steps 1, 1a, 1b
    vpt       *_micron_space.parquet                          VPT or Vizgen cellpose, step 1b
    cellpose  labels.npy or *-full-merged-labels-compressed.npz   the lab's merged cellpose, step 1a

--steps narrows the files counted as needed to the steps you mean to run; which ones a file serves
is printed beside it. --raw also keeps the rest of what the instrument and the segmentation run
wrote -- every z plane, the .vzg, the per-FOV cellpose labels -- but not what was derived later. --globus writes those files as a batch for `globus transfer --batch`, each
path relative to --root on both sides, so the tree keeps its shape on the far end.
"""

import argparse
import csv
import fnmatch
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path, PurePosixPath

# The plane merscope() reads images from -- Z_LAYER in create_spatialdata*.py.
Z_LAYER = 3

# kind -> [(label, patterns, steps)]. Patterns are alternatives in order of preference: the first
# that matches anything is the one taken, as step 1 takes the boundary parquet over the per-FOV
# HDF5 and step 1a the compressed labels over labels.npy. They are globs relative to the directory, the glob in the last component only, and one ending in "/"
# matches whole directories, sent recursively. Step 1b takes only the images and transform from a
# region: its counts, metadata and boundaries come from vpt_outputs, and its transcripts are staged
# over the region's own.
#
# "raw" is not a step: it is what --raw adds, the rest of what the instrument or the segmentation
# run wrote, so a transfer keeps the source data and not only the pipeline's inputs. What the lab
# derived afterwards -- z stacks, 8-bit composites, transcript assignments, Baysor runs -- is not
# raw and is left out.
RULES = {
    "region": [
        ("counts", ["cell_by_gene.csv"], {"1", "1a"}),
        ("metadata", ["cell_metadata.csv"], {"1"}),
        ("boundaries", ["cell_boundaries.parquet", "cell_boundaries/feature_data_*.hdf5"], {"1"}),
        ("transcripts", ["detected_transcripts.csv"], {"1", "1a"}),
        ("transform", ["images/micron_to_mosaic_pixel_transform.csv"], {"1", "1a", "1b"}),
        ("mosaic images", [f"images/mosaic_*_z{Z_LAYER}.tif"], {"1", "1a", "1b"}),
        ("other z planes", ["images/mosaic_*_z*.tif"], {"raw"}),
        ("image manifest", ["images/manifest.json"], {"raw"}),
        ("vzg", ["*.vzg"], {"raw"}),
        ("summary", ["summary.png"], {"raw"}),
    ],
    "vpt": [
        ("counts", ["cellpose_cell_by_gene.csv"], {"1b"}),
        ("metadata", ["cellpose_cell_metadata.csv"], {"1b"}),
        ("boundaries", ["cellpose_micron_space.parquet"], {"1b"}),
        ("transcripts", ["detected_transcripts.csv"], {"1b"}),
    ],
    "cellpose": [
        # The same array, ~1% the size; step 1a streams it a plane at a time.
        ("labels", ["*-full-merged-labels-compressed.npz", "labels.npy"], {"1a"}),
        ("label map", ["*-cellp-label-map.npz"], {"1a"}),
        ("counts", ["*-cell-by-gene.npz"], {"1a"}),
        ("bounding boxes", ["*-slicee.pkl"], {"1a"}),
        ("centres", ["cell-coms.npz"], {"1a"}),
        ("volumes", ["cell-vols.npz"], {"1a"}),
        ("per-FOV labels", ["[0-9]*/"], {"raw"}),
        ("FOV maps", ["fov_dict*.pkl"], {"raw"}),
        ("label maps", ["map_dict.pkl"], {"raw"}),
        ("parameters", ["cellpose_params.json"], {"raw"}),
        ("region used", ["region_dir.txt"], {"raw"}),
        ("logs", ["logs/"], {"raw"}),
    ],
}
COLUMN_KIND = {"path": "region", "vpt_path": "vpt", "cellpose_path": "cellpose"}
ALL_STEPS = {"1", "1a", "1b"}
# Image and per-FOV directories hold thousands of files and never another dataset.
PRUNE = {"images", "cell_boundaries", "result_tiles"}
COLLECTION = re.compile(r"^([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}):(.*)$")


class Local:
    name = "local"

    def __init__(self):
        self.cache = {}

    def listdir(self, path):
        """Return [(name, is_dir, size)], or None when the directory cannot be read."""
        if path not in self.cache:
            try:
                entries = []
                for e in os.scandir(path):
                    try:
                        is_dir = e.is_dir()  # follows links, as the pipeline does
                        entries.append((e.name, is_dir, 0 if is_dir else e.stat().st_size))
                    except OSError:
                        pass  # a broken link
                self.cache[path] = entries
            except OSError as err:
                print(f"    !! {err}")
                self.cache[path] = None
        return self.cache[path]

    def du(self, path):
        total = 0
        for dirpath, _, filenames in os.walk(path, followlinks=True):
            for f in filenames:
                try:
                    total += os.stat(os.path.join(dirpath, f)).st_size
                except OSError:
                    pass
        return total

    def key(self, path):
        return os.path.realpath(path)


class Globus:
    """A collection, listed one directory per `globus ls` call."""

    def __init__(self, collection):
        self.collection = collection
        self.name = collection
        self.cache = {}
        self.cli = shutil.which("globus")
        if not self.cli:
            sys.exit("A collection path needs the globus CLI on PATH -- conda activate globus.")

    def listdir(self, path):
        if path not in self.cache:
            result = subprocess.run(
                [self.cli, "ls", "-F", "json", f"{self.collection}:{str(path).rstrip('/')}/"],
                capture_output=True, text=True,
            )
            if result.returncode:
                print(f"    !! globus ls {path}: {result.stderr.strip()}")
                self.cache[path] = None
            else:
                self.cache[path] = [
                    (e["name"], e["type"] == "dir", 0 if e["type"] == "dir" else e["size"])
                    for e in json.loads(result.stdout)["DATA"]
                ]
        return self.cache[path]

    def du(self, path):
        return None  # a recursive listing would cost one call per directory

    def key(self, path):
        return str(path)


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n} B" if unit == "B" else f"{n:,.1f} {unit}"
        n /= 1024


def classify(fs, directory, entries):
    names = {n for n, _, _ in entries}
    dirs = {n for n, is_dir, _ in entries if is_dir}
    if "images" in dirs:
        images = fs.listdir(directory / "images") or []
        if any(n == "micron_to_mosaic_pixel_transform.csv" for n, _, _ in images):
            return "region"
    if "labels.npy" in names or any(n.endswith("-full-merged-labels-compressed.npz") for n in names):
        return "cellpose"
    if any(n.endswith("_micron_space.parquet") for n in names):
        return "vpt"
    return None


def discover(fs, root):
    """Yield (kind, directory) for everything under root the pipeline can read."""
    stack = [root]
    while stack:
        directory = stack.pop()
        entries = fs.listdir(directory)
        if entries is None:
            continue
        kind = classify(fs, directory, entries)
        if kind:
            yield kind, directory
        if kind == "cellpose":
            continue  # hundreds of per-FOV directories, none of them a dataset
        subdirs = sorted((n for n, is_dir, _ in entries if is_dir and n not in PRUNE), reverse=True)
        stack.extend(directory / n for n in subdirs)


def from_samplesheet(sheet):
    with open(sheet, newline="") as fh:
        for row in csv.DictReader(fh):
            for column, kind in COLUMN_KIND.items():
                if row.get(column):
                    yield kind, Path(row[column]).absolute(), row["sample"]


def unread_size(fs, path, matched):
    """(bytes, directories not sized) under path outside matched, descending only where matched."""
    if not any(path in m.parents for m in matched):
        n = fs.du(path)
        return (0, 1) if n is None else (n, 0)
    total, unsized = 0, 0
    for name, is_dir, n in fs.listdir(path) or []:
        child = path / name
        if is_dir:
            n, u = unread_size(fs, child, matched)
            total, unsized = total + n, unsized + u
        elif child not in matched:
            total += n
    return total, unsized


def inventory(fs, kind, directory, steps):
    """Print one directory's report; return (needed [(path, size, is_dir)], missing count).

    A directory's size is None when the filesystem cannot say without walking it.
    """
    top = fs.listdir(directory)
    if top is None:
        print("    !! not a readable directory (missing, or a broken link)")
        return [], 1

    needed, matched, missing, read_bytes = [], set(), 0, 0
    for label, patterns, used_by in RULES[kind]:
        hits, dirs = {}, set()
        subdirs = {n for n, is_dir, _ in top if is_dir}
        for pat in patterns:
            if pat.endswith("/"):
                for name in sorted(subdirs):
                    if fnmatch.fnmatchcase(name, pat[:-1]) and directory / name not in matched:
                        hits[directory / name] = fs.du(directory / name)
                        dirs.add(directory / name)
                continue
            parent, _, glob = pat.rpartition("/")
            if parent and parent not in subdirs:
                continue  # an alternative form this directory does not use
            where = directory / parent if parent else directory
            for name, is_dir, n in (fs.listdir(where) if parent else top) or []:
                path = where / name
                if not is_dir and fnmatch.fnmatchcase(name, glob) and path not in matched:
                    hits[path] = n
            if hits:
                break
        raw = used_by == {"raw"}
        tag = "raw" if raw else "step " + ", ".join(sorted(used_by))
        wanted = bool(used_by & steps)
        if not hits:
            missing += wanted and not raw
            status = "MISSING" if wanted and not raw else "absent"
            print(f"    {status:<8} {label:<15} {' or '.join(patterns):<50} {'':>10}  {tag}")
            continue
        sizes = [n for n in hits.values() if n is not None]
        total = sum(sizes)
        if not raw:
            read_bytes += total
        matched.update(hits)
        if wanted:
            needed.extend((path, n, path in dirs) for path, n in hits.items())
        noun = "directories" if dirs else "files"
        shown = next(iter(hits)).relative_to(directory) if len(hits) == 1 else f"{len(hits)} {noun}"
        size = human(total) if len(sizes) == len(hits) else "not sized" if not sizes else f"{human(total)}+"
        print(f"    {'need' if wanted else 'skip':<8} {label:<15} {str(shown):<50} {size:>10}  {tag}")

    unread, unsized = [], 0
    for name, is_dir, n in top:
        path = directory / name
        if path in matched:
            continue
        if is_dir:
            n, u = unread_size(fs, path, matched)
            unsized += u
            if n == 0 and u == 0:
                continue
        unread.append((name + ("/" if is_dir else ""), n))
    if unread:
        extra = f", plus {unsized} director{'y' if unsized == 1 else 'ies'} not sized" if unsized else ""
        print(f"    left out, {human(sum(n for _, n in unread))}{extra}:")
        unread.sort(key=lambda e: -e[1])
        for name, n in unread[:8]:
            print(f"      {human(n):>10}  {name}")
        if len(unread) > 8:
            print(f"      {human(sum(n for _, n in unread[8:])):>10}  ... {len(unread) - 8} more entries")
    print(f"    needed: {human(sum(n or 0 for _, n, _ in needed))} ({human(read_bytes)} the pipeline can read)")
    return needed, missing


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs", nargs="+",
                        help="directories to walk, COLLECTION_ID:/path to walk over globus, or samplesheets")
    parser.add_argument("--steps", default="1,1a,1b",
                        help="comma-separated steps to keep files for (default: %(default)s)")
    parser.add_argument("--raw", action="store_true",
                        help="also keep the rest of the instrument's and the segmentation run's own output")
    parser.add_argument("--globus", type=Path, help="write the needed files here as a globus --batch file")
    parser.add_argument("--root", help="base path the batch lines are relative to (default: the needed files' common path)")
    args = parser.parse_args()

    steps = set(args.steps.split(","))
    if not steps <= ALL_STEPS:
        parser.error(f"--steps takes {', '.join(sorted(ALL_STEPS))}, got {args.steps}")
    if args.raw:
        steps.add("raw")

    local, remotes, found = Local(), {}, []
    for item in args.inputs:
        match = COLLECTION.match(item)
        if match:
            fs = remotes.setdefault(match[1], Globus(match[1]))
            found.extend((fs, kind, d, None) for kind, d in discover(fs, PurePosixPath(match[2] or "/")))
        elif item.endswith(".csv"):
            found.extend((local, kind, d, sample) for kind, d, sample in from_samplesheet(item))
        else:
            found.extend((local, kind, d, None) for kind, d in discover(local, Path(item).absolute()))

    needed, seen, missing, by_kind = {}, set(), 0, {}
    for fs, kind, directory, sample in found:
        key = (fs.name, kind, fs.key(directory))
        if key in seen:
            continue
        seen.add(key)
        where = "" if fs is local else f"{fs.name}:"
        print(f"\n[{kind}] {where}{directory}" + (f"  (sample {sample})" if sample else ""))
        files, m = inventory(fs, kind, directory, steps)
        missing += m
        count, n = by_kind.get(kind, (0, 0))
        by_kind[kind] = (count + 1, n + sum(s or 0 for _, s, _ in files))
        # A file reached by two routes -- a region named in two samplesheets -- is sent once.
        for path, size, is_dir in files:
            needed[(fs.name, path)] = (size, is_dir)

    total = sum(size or 0 for size, _ in needed.values())
    unsized = sum(size is None for size, _ in needed.values())
    print(f"\n== Summary for steps {', '.join(sorted(steps))}")
    for kind, (count, n) in sorted(by_kind.items()):
        print(f"  {count:>3} {kind:<9} directories  {human(n):>10} needed")
    ndirs = sum(is_dir for _, is_dir in needed.values())
    print(f"  {len(needed) - ndirs:,} files and {ndirs:,} whole directories, {human(total)} to transfer"
          + (f", plus {unsized:,} directories not sized" if unsized else ""))
    if missing:
        print(f"  {missing} needed file(s) missing -- see MISSING above.")

    links = [p for (where, p) in needed if where == "local" and os.path.islink(p)]
    if links:
        print(f"  {len(links)} needed file(s) are symlinks; globus does not follow them. First: {links[0]}")

    if args.globus and needed:
        sources = {where for where, _ in needed}
        if len(sources) > 1:
            sys.exit("  --globus: the needed files span more than one source; run once per source.")
        (source,) = sources
        paths = sorted(p for _, p in needed)
        root = PurePosixPath(args.root or os.path.commonpath([str(p) for p in paths]))
        with open(args.globus, "w") as fh:
            for p in paths:
                rel = shlex.quote(str(PurePosixPath(p).relative_to(root)))
                recursive = "--recursive " if needed[(source, p)][1] else ""
                fh.write(f"{recursive}{rel} {rel}\n")
        src = f"{source}:{root}" if source != "local" else f"SRC_COLLECTION_ID:{root}"
        print(f"\n  wrote {args.globus}, paths relative to {root}. Transfer with:")
        print(f"    globus transfer {src} DST_COLLECTION_ID:/destination/path --batch {args.globus} --label merfish-minimal")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
