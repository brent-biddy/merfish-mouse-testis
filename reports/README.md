# Reports

One directory per render. Each holds the document, its figures and its deck together, and
is self-contained — the links inside a document are relative to its own directory, so a
directory can be renamed or moved without breaking anything in it.

The markdown and the figures are committed, being what GitHub renders. The deck is not,
being the same content in a format git cannot diff; rebuild it by re-rendering.

A render publishes into `reports/<qmd>_<run_id>_<to>/`, one format per render, and a
`render_sample` run nests one directory per sample inside that. `--to` defaults to `pptx`,
which is gitignored; `--to gfm` produces the copy this directory commits.

A render with a fresh run id cannot displace an earlier one, the timestamp being unique per
invocation — so this directory accumulates: prune the renders not worth keeping before
committing them. Passing `--run_id` gives that up: two runs sharing a run id and format
publish into one directory, and `publishDir` copies over what is there.

## Index

Written by hand. A render cannot say what it was for, and a pipeline that wrote this file
would overwrite the part that matters. Add a row when you commit a report.

| Report | Samples | Segmentation | What it was for |
|--------|---------|--------------|-----------------|
| [celltype_report_b2r0_report](celltype_report_b2r0_report_gfm/b2r0_cellpose3d/celltype_report.md) | b2r0_cellpose3d | cellpose (lab's own 2023 pipeline, see `b2r0_cellpose3d` provenance) | First celltype report on a real sample, not the `u2os_test` fixture |

<!-- Row format — the link points at the document inside the directory:
| [celltype_report_cellpose_cmp](celltype_report_cellpose_cmp/celltype_report.md) | testis_01 | cellpose | Comparison against the Vizgen boundaries before settling on one |
-->
