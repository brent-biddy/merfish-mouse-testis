# Project notes

Where the 16-section run stands. Dated 2026-10-30; update it as jobs finish. On OSCER the repo
is `/ourdisk/hpc/lilab/babiddy/dont_archive/merfish-mouse-testis`, and the data is beside it in
`merfish-mouse-testis-data` and under `alex_vizgen/ahs_vizgen`.

## Goal

Cluster every tissue section, then call cell types with the centroid correlation in step 3.
Alex Vargo's dataset has 16 sections; the pipeline has run 14 of them.

## Status by section

| Sections | Segmentation | Step 1a | Step 1b | Step 2 cluster | Step 3 annotate |
|---|---|---|---|---|---|
| 2ar0 2ar1 2br0 2br1 a1r1 b2r0 b2r1 b2r2 c1r0 c1r2 d2r0 d2r1 d2r2 | 3D cellpose, in `merfish-mouse-testis-data/segmentation/` | done | done | done | not run |
| c1r1 | same | done | queued | queued, chained after 1b | not run |
| a1r0-upper, a1r0-lower | in `alex_vizgen/ahs_vizgen/cellposetest/`, one raster for both | split queued | not started | not started | not run |

Outputs for the 13 finished sections are in `results/cellpose3d_all/`. Its handoff sheets list
c1r1 too, but its store and table were never written, so they cannot be reused for it.

## Jobs submitted 2026-10-30, both pending on Priority when last checked

- **c1r1 steps 1b and 2**: Slurm job 35614280, run id `c1r1_finish`. A detached Nextflow
  process chains step 2 after step 1b; if it dies, re-run with `-resume --run_id c1r1_finish`.
  Log `c1r1_finish.out`, results in `results/c1r1_finish/`. Sheet: `results/cellpose3d_all/c1r1_create.csv`.
- **a1r0 split**: Slurm job 35614179, log `a1r0_split_35614179.log`, output
  `merfish-mouse-testis-data/segmentation/a1r0_split/`. It runs `bin/split_cellpose_section.py`.

## a1r0

Two tissue sections sit in one region, so one merge.py raster and one counts matrix cover both.
Facts, from Alex's `section.py` and `xscript-metafile-creation.py`:

- The cut is `global_y` 3520 microns, which is mosaic row 33,500. Lower is `global_y > 3520`;
  upper is the rest. Cells split the same way: 84,563 above, 157,060 below, 241,623 in all.
- The raster is cropped. Its label offset is `-245,-290`; prep's own offset check passes for all
  cells, and `0,0` fails. The lower section's raster starts at the cut, so its offset is `33500,-290`.
- `split_cellpose_section.py` writes a segmentation directory and a region copy per section, with
  only that section's transcripts, so step 1a needs no change.

When the split job finishes: check its printed counts against the numbers above, then add rows
to `assets/samplesheet_cellpose_prep.csv`, one per section, with `path` set to
`<outdir>/regions/<name>`, `cellpose_path` to `<outdir>/<name>`, and the offset above. Then run
steps 1a, 1b and 2 on them. Neither the split script nor the filter on transcripts has been run yet.

## Annotation references

`annotate_celltypes.py --reference` takes either table in `assets/reference/`:

- `shami_human_testis_centroids.csv.gz`, the default: human, 11 major types, ln(mean + 1).
- `green_mouse_testis_centroids.csv.gz`: mouse, 32 types (germ stages, Sertoli subtypes), the 491
  panel genes. Alex's `VizgenFinal491_centroids_rawRename7.txt`. His `greenMeta.py` attributes it
  to Green et al. 2018. **Open question for Alex:** how its values were computed. His own script
  uses Pearson on log2 of counts per 1000, ours is Spearman, so whether it matches our centroid
  recipe is unconfirmed.

## Open items

- Step 3 has not been run on any section; choose the reference.
- The OSCER checkout's `origin` is `merfish_testis.git`, not `merfish-mouse-testis.git`. Pulls
  work, but it is probably stale.
- `CLAUDE.md` says c1r1 is the only cropped raster of the fourteen; a1r0 is a second.
- Nothing from Alex's UMich storage is needed so far; everything used is already on OSCER.

## Working on OSCER

Open a shared connection from your own terminal: `ssh-add ~/.ssh/id_ed25519` then `ssh -fN oscer`.
It closes after 30 idle minutes; `ssh -O exit oscer` closes it sooner. Keep heavy work off the
login node and run it through Slurm. Nextflow comes from `module load Nextflow/24.10.2`
(or your conda env).
