"""Rebuild ``data/pbmcs3k_pre.h5ad``, the PBMC3k matrix the celltype arm reads.

``celltype/config.yaml`` points at ``data/pbmcs3k_pre.h5ad`` and nothing in the
repository built it. The nearest thing was cell 4 of
``clustering/legacy/resolution_sweep.ipynb``, which writes a file of that name --
but not this one: it writes raw counts, before normalisation, with a single
``counts`` layer. The committed artefact holds log1p-normalised values in ``X``,
carries a second ``log1p_norm`` layer, and its ``.raw`` is log1p over all 13 714
genes. So the notebook records an earlier version of the pipeline.

The order below was recovered from the artefact itself rather than from the
notebook, and each step is pinned by something measurable in it:

* ``layers['counts']`` is integral, so counts are captured *before* normalising;
* ``log1p_norm`` equals ``log1p(counts / n_counts * 1e4)`` where ``n_counts`` is
  the **full-gene** total, not the HVG-subset total -- so ``normalize_total``
  ran while all 13 714 genes were still present;
* ``.raw`` is non-integral and has 13 714 genes, so it was assigned *after*
  log1p and *before* the HVG subset;
* ``X`` equals the ``log1p_norm`` layer exactly.

Run from ``reproducibility/``::

    python -m celltype.preprocess            # writes data/pbmcs3k_pre.h5ad
    python -m celltype.preprocess --check    # rebuild and diff against the existing file
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import scanpy as sc

from paths import require_input

# The PBMC3k matrix is a public 10x download; clustering/config.yaml points at the
# tarball and clustering/legacy/resolution_sweep.ipynb records the URL.
TENX_URL = (
    "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/"
    "pbmc3k_filtered_gene_bc_matrices.tar.gz"
)
SHARED_DATA = pathlib.Path(__file__).resolve().parent.parent / "data"
MTX_DIR = SHARED_DATA / "filtered_gene_bc_matrices" / "hg19"
OUT = SHARED_DATA / "pbmcs3k_pre.h5ad"

N_TOP_GENES = 2000
TARGET_SUM = 1e4


def build():
    mtx = require_input(
        MTX_DIR,
        what="the raw PBMC3k 10x matrix (barcodes.tsv, genes.tsv, matrix.mtx)",
        source=(
            f"public 10x download. Fetch and unpack with:\n"
            f"           curl -o data/pbmc3k.tar.gz {TENX_URL}\n"
            f"           tar -xzf data/pbmc3k.tar.gz -C data/"
        ),
    )
    adata = sc.read_10x_mtx(str(mtx), var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()

    # --- QC, the standard PBMC3k tutorial thresholds -----------------------
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    mito = adata.var_names.str.startswith("MT-")
    counts_per_cell = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["percent_mito"] = (
        np.asarray(adata[:, mito].X.sum(axis=1)).ravel() / counts_per_cell
    )
    adata.obs["n_counts"] = counts_per_cell
    adata = adata[adata.obs["n_genes"] < 2500, :]
    adata = adata[adata.obs["percent_mito"] < 0.05, :].copy()

    # --- counts kept before anything rescales them -------------------------
    adata.layers["counts"] = adata.X.copy()

    # --- HVG on the counts, while all genes are still present --------------
    # seurat_v3 expects raw counts, which is why it reads the layer.
    sc.pp.highly_variable_genes(
        adata, n_top_genes=N_TOP_GENES, layer="counts", flavor="seurat_v3"
    )

    # --- normalise over the full gene set, then log1p ----------------------
    # Doing this before the HVG subset is what makes log1p_norm's denominator
    # the full-gene total; subsetting first would divide by the HVG total and
    # give a different matrix.
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)

    # --- raw: log1p over all genes -----------------------------------------
    adata.raw = adata

    # --- subset to the highly variable genes -------------------------------
    adata = adata[:, adata.var["highly_variable"]].copy()
    adata.layers["log1p_norm"] = adata.X.copy()
    return adata


def _dense(m):
    return m.toarray() if hasattr(m, "toarray") else np.asarray(m)


def compare(built, existing_path):
    """Report how the rebuild differs from the committed artefact."""
    old = sc.read_h5ad(existing_path)
    ok = True

    def check(label, cond, detail=""):
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'ok  ' if cond else 'DIFF'} {label}{detail}")

    check("shape", built.shape == old.shape, f"  {built.shape} vs {old.shape}")
    check("obs columns", list(built.obs.columns) == list(old.obs.columns))
    check("var names", list(built.var_names) == list(old.var_names))
    check("obs names", list(built.obs_names) == list(old.obs_names))
    check("layers present", set(built.layers) == set(old.layers))
    for key in ("counts", "log1p_norm"):
        if key in built.layers and key in old.layers:
            a, b = _dense(built.layers[key]), _dense(old.layers[key])
            m = np.abs(a - b).max()
            check(f"layer {key}", m < 1e-5, f"  max|diff| {m:.3g}")
    a, b = _dense(built.X), _dense(old.X)
    m = np.abs(a - b).max()
    check("X", m < 1e-5, f"  max|diff| {m:.3g}")
    if built.raw is not None and old.raw is not None:
        check("raw shape", built.raw.shape == old.raw.shape,
              f"  {built.raw.shape} vs {old.raw.shape}")
        m = np.abs(_dense(built.raw.X) - _dense(old.raw.X)).max()
        check("raw X", m < 1e-5, f"  max|diff| {m:.3g}")
    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--check", action="store_true",
                   help="rebuild and compare against the existing file; write nothing")
    p.add_argument("--out", default=str(OUT), help=f"output path (default: {OUT})")
    args = p.parse_args()

    adata = build()
    print(f"Built {adata.shape[0]} cells x {adata.shape[1]} genes "
          f"(raw {adata.raw.shape[1]} genes)")

    if args.check:
        target = pathlib.Path(args.out)
        if not target.exists():
            raise SystemExit(f"--check needs an existing file to compare against: {target}")
        print(f"Comparing against {target}:")
        sys.exit(0 if compare(adata, target) else 1)

    adata.write_h5ad(args.out)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
