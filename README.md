# LN's $t$-test

Differential expression testing on an asymptotically unbiased log-fold-change estimator.

The repository holds two things:

- **`src/lntest/`** — the method, as an installable package. `ln_test.py` is the estimator; `scanpy_wrapper.py` exposes it as `rank_genes_groups_ln`, shaped like `scanpy.tl.rank_genes_groups`.
- **`reproducibility/`** — every result in the Nature submission that depends on code, one directory per experiment.

## Install

```bash
pip install -e .
```

That is all that is needed. `lntest` depends on `numpy`, `scipy` and `statsmodels` — the last for Benjamini-Hochberg, imported lazily and only when you ask for it, which is the arrangement scanpy itself uses.

It does **not** depend on `anndata` or `scanpy`. `rank_genes_groups_ln` is duck-typed on the object you hand it, so it works with an `AnnData` without the package needing to import one.

## Use

```python
import scanpy as sc
from lntest import rank_genes_groups_ln

rank_genes_groups_ln(adata, groupby="leiden", layer="norm_counts")
adata.uns["rank_genes_groups"]["logfoldchanges"]
adata.uns["rank_genes_groups"]["lfc_se"]      # SE of the LFC; lfc ± 1.96·se is the interval
```

Takes **normalised, not log-transformed** counts. The interval is the method's distinguishing claim, so `lfc_se` is returned alongside the estimate.

Two options are worth knowing about:

- `corr_method=` — `"benjamini-hochberg"` or `"bonferroni"`, scanpy's own vocabulary and scanpy's
  implementation. The default is `"bonferroni"` rather than scanpy's `"benjamini-hochberg"`, because
  this function was hardcoded to Bonferroni before the parameter existed and every published number
  was produced that way.
- `trigamma=` — `"exact"` (ψ₁, the default) or `"recomb25"`, the `1/x` approximation every published RECOMB number was produced with. Which is correct is an open question; the reproducibility scripts request `"recomb25"` explicitly so that tree reproduces the submission by construction.

## Reproducing the Nature submission

**See [`reproducibility/README.md`](reproducibility/README.md).** It is the only reproduction instruction that is maintained: the environment, the run convention, what each arm needs, and what is known not to work.

Two things that catch people immediately:

- run everything from `reproducibility/` as `python -m <arm>.<script>`, not by file path;
- `export R_HOME="$CONDA_PREFIX/lib/R"`, or three of the arms abort in rpy2 before running any of your code.

## Tests

```bash
pytest
```

Covers the package only. **`reproducibility/` has no test suite** — nothing there checks a change.

## Contributing

`main` here is a shared branch and this repository is a fork; work happens on `*_dev` branches and reaches `main` by pull request. Method choices that someone will later have to justify are recorded as decision records rather than left in commit messages.
