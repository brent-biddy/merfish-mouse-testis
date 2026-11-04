#!/usr/bin/env python3
"""
build_green_reference.py - Build the Green et al. 2018 reference tables over a gene panel.

Writes the two reference tables the per-section notebook loads and never builds, each
genes x types over the panel genes only:

    green_major.csv   the 11 major types, ln(mean TP10K + 1), rebuilt from the published
                      counts and per-cell labels
    green_gc.csv      the 12 germ-cell clusters GC1-GC12, as published

The 11 major types are rebuilt because GEO has the counts and labels but not the centroids:
every cell's counts are scaled by its total over all genes, to TP10K, then averaged per type
and log1p'd. The DGE is 37,241 genes x 34,633 cells, so it is read in chunks and takes a few
minutes. The GC1-12 centroids are already log(average normalized expression + 1) and are used
as they are.

Panel genes renamed since the reference was made are looked up under the name the reference
uses, so they are not dropped. A panel gene found under neither name is left out and reported.

Inputs, from GEO GSE112393, in --geo_dir:

    GSE112393_MergedAdultMouseST25_DGE.txt.gz
    GSE112393_MergedAdultMouseST25_PerCellAttributes.txt.gz
    GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx

Usage:
    build_green_reference.py \\
        --geo_dir data/reference/green2018 \\
        --panel_table data/processed/b2r0/b2r0_cellpose3d.create_spatialdata_cellpose.h5ad \\
        --outdir results/section_analysis/reference
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# Per-cell label codes in the PerCellAttributes file.
CELL_TYPE_CODES = {
    1: "InnateLymph", 2: "Macrophage", 3: "Endothelial", 4: "Myoid", 5: "Leydig",
    6: "Sertoli", 7: "Unknown", 8: "SPG", 9: "Scytes", 10: "STids", 11: "Elongating",
}

# Panel genes the reference knows under an older name: panel name -> name in the reference.
RENAMES = {
    "Plpp3": "Ppap2b", "Meioc": "1700012B09Rik", "Mesd": "Mesdc2",
    "Prxl2a": "Fam213a", "Ints6l": "Ddx26b", "Map7d3": "Mtap7d3",
}

DGE = "GSE112393_MergedAdultMouseST25_DGE.txt.gz"
ATTRIBUTES = "GSE112393_MergedAdultMouseST25_PerCellAttributes.txt.gz"
GERM_CLUSTERS = "GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx"


def parse_args():
    parser = argparse.ArgumentParser(description="Build the Green 2018 reference tables over a gene panel")
    parser.add_argument("--geo_dir", required=True, help="Directory holding the three GEO files")
    parser.add_argument("--panel_table", required=True, help="An .h5ad whose var_names are the panel genes")
    parser.add_argument("--outdir", default=".", help="Directory to write the two tables into")
    parser.add_argument("--chunksize", type=int, default=2000, help="DGE rows read at a time")
    return parser.parse_args()


def match_panel_genes(panel):
    """Panel gene -> the lower-cased name to look for in the reference, renamed genes included."""
    return {gene: RENAMES.get(gene, gene).lower() for gene in panel}


def report_left_out(panel, kept, label):
    left_out = [gene for gene in panel if gene not in kept]
    print(f"{label}: {len(kept)} of {len(panel)} panel genes; left out: {left_out}")


def build_major(panel, geo_dir, chunksize):
    """Mean log TP10K per type over the panel genes, from the GEO counts and per-cell labels."""
    # '#' starts the header line, so the labels begin after four lines of comments.
    labels = pd.read_csv(
        geo_dir / ATTRIBUTES, sep="\t", skiprows=4, header=None,
        usecols=[0, 2], names=["cell", "code"], index_col="cell",
    )
    wanted = match_panel_genes(panel)
    wanted_lower = set(wanted.values())

    cell_total, rows = None, []
    for chunk in pd.read_csv(geo_dir / DGE, sep="\t", index_col=0, chunksize=chunksize):
        totals = chunk.sum(axis=0)  # every gene counts towards a cell's total
        cell_total = totals if cell_total is None else cell_total + totals
        rows.append(chunk[chunk.index.astype(str).str.lower().isin(wanted_lower)])

    counts = pd.concat(rows)
    counts.index = counts.index.astype(str).str.lower()
    counts = counts.loc[~counts.index.duplicated()]

    # Back to the panel's own spelling, keeping the genes the reference has.
    found = {gene: name for gene, name in wanted.items() if name in counts.index}
    counts = counts.loc[list(found.values())].set_axis(list(found.keys()))
    report_left_out(panel, counts.index, "green_major")

    tp10k = counts.div(cell_total.loc[counts.columns], axis=1) * 1e4
    types = labels["code"].map(CELL_TYPE_CODES).loc[tp10k.columns]
    return np.log1p(tp10k.T.groupby(types.to_numpy()).mean().T)  # genes x types


def build_germ_clusters(panel, geo_dir):
    """GC1-GC12 over the panel genes, as published."""
    # The sheet opens with description lines; the header is row 5, then nGene and nCell rows.
    sheet = pd.read_excel(geo_dir / GERM_CLUSTERS, sheet_name="Average", header=5, index_col=0)
    sheet.columns = [str(c) for c in sheet.columns]
    sheet = sheet.drop(["nGene", "nCell"]).astype(float)

    by_lower = {name.lower(): name for name in sheet.index.astype(str)}
    found = {gene: by_lower[name] for gene, name in match_panel_genes(panel).items() if name in by_lower}
    table = sheet.loc[list(found.values())].set_axis(list(found.keys()))
    report_left_out(panel, table.index, "green_gc")
    return table


def main():
    args = parse_args()
    geo_dir, outdir = Path(args.geo_dir), Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    panel = list(ad.read_h5ad(args.panel_table, backed="r").var_names)

    build_major(panel, geo_dir, args.chunksize).to_csv(outdir / "green_major.csv")
    build_germ_clusters(panel, geo_dir).to_csv(outdir / "green_gc.csv")
    print(f"wrote {outdir / 'green_major.csv'} and {outdir / 'green_gc.csv'}")


if __name__ == "__main__":
    main()
