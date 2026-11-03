# Notebook design

How the per-section Quarto notebook is built. `PLAN.md` is what the analysis does and in what
order; this file is how the notebook carries it out. Only Phases B and C are designed so far.

## One object

**One object passes through the whole notebook.** Section 1 loads the section's table as a
single AnnData, and every later section adds to that same object rather than making a new one:
raw counts kept in `layers["counts"]` before normalization overwrites `X`, embeddings in
`obsm`, clusters, scores and calls in `obs`, markers and parameters in `uns`. A checkpoint is
that object written to disk, and loading one puts the notebook back where it was. The
SpatialData store's images and boundaries are not loaded into it; tissue figures read the
shapes they need from the store when they draw.

Memory clean-up is not planned. Adding to one object in place leaves no stale copies to free,
and the object is small: about 147,000 cells by 491 genes is under 300 MB dense, plus the
neighbour graph. Two habits keep it that way: no `adata.copy()` per section, and anything big
that is only needed briefly, such as a full reference dataset, is deleted once used.

## Checkpoints

Each section caches its result as an on-disk checkpoint: a file (`.h5ad`, `.parquet`) named
by a key built from that section's parameters and the key of the section before it. A section
whose file exists loads it instead of recomputing, so editing a later section never reruns an
earlier one, and changing a parameter reruns only that section and the ones after it. A code
change the parameters do not capture needs the checkpoint deleted, or the section's code
hashed into its key. The checkpoints are ordinary files, readable outside the notebook.

The checkpoints are the notebook's intermediate results, written into `results/` like the
pipeline's and gitignored with it:

```
results/<notebook>/<sample>/<section>/<sample>.<section>.<key>.<ext>
```

The layout is the pipeline's `results/<run_id>/<sample>/<step>/` with the notebook's name in
place of the run id. A timestamped run id would start every render in an empty directory and
find nothing to load; a fixed one lets each render reuse what the last one wrote. The key in
the filename keeps results from different parameters side by side rather than overwriting
each other.

## Sections

### Phase B, clustering

1. Load the section's table
2. Filter cells
3. Normalize
4. PCA, neighbours and UMAP
5. Leiden clustering at 0.5, 1.0, 1.5 and 2.0

   **Checkpoint 1, clustered:** the filtered table with a cluster column per resolution.

Changing anything here reruns sections 2 to 5, and invalidates Phase C's checkpoint 2 with it,
since that key includes checkpoint 1's. Marker genes are not computed here: they depend on
the resolution being annotated, so they are Phase C's first section.

A third checkpoint could go after normalization (section 3), so that changing PCA, neighbours
or Leiden does not redo filtering and normalization. Not needed yet; add it if clustering
b2r0 turns out slow.

### Phase C, annotation

6. Choose the resolution to annotate, 1.0 to start
7. Marker genes per cluster, at that resolution
8. Per-cell correlation against the chosen reference
9. Cluster composition of the per-cell calls
10. Cluster calls

    **Checkpoint 2, annotated:** the whole table, with everything sections 1 to 10 added: the
    clusters at every resolution, the markers, every `corr_<type>` score, the per-cell calls and
    the cluster calls. The later phases read whichever columns they need from it, so which
    label each one uses is decided when that phase is built.

Changing the resolution, the markers, the correlation or the calls reruns only sections 6 to
10, starting from checkpoint 1, so annotation is iterated on without reclustering. The
resolution, the reference and the purity threshold are all part of checkpoint 2's key, so a run
with each reference leaves its own checkpoint 2 side by side.

## CPU or GPU

**The notebook starts on CPU (scanpy).** b2r0 has about 147,000 cells after filtering, over
491 genes, which CPU clusters in minutes, and a CPU job runs on any OSCER node or locally with
little or no wait.

GPU (rapids-singlecell, as step 2 used) is the option if clustering becomes the slow part. It
clusters faster, but OSCER's GPU queue is the wait: step 2's jobs on `sooner_gpu_test` were
scheduled hours apart, the last a day out. The `wsl` machine's own GPU has no queue.
rapids-singlecell mirrors scanpy's API, so moving over is a switch in the notebook rather than
a rewrite. The two Leiden implementations do not give identical clusters, so the backend is
part of checkpoint 1's key.
