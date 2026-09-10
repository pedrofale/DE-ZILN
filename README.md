# LN's $t$-test

Differential expression testing on an asymptotically unbiased log-fold-change estimator.

The repository holds two things:

- **`src/lntest/`** — the method, as an installable package. `rank_genes_groups_ln` is shaped like `scanpy.tl.rank_genes_groups`; `get_LN_lfcs` is the estimator underneath. Both import from `lntest` directly — the implementation modules are private.
- **`reproducibility/`** — every result in the paper that depends on code, one directory per experiment. See [`reproducibility/README.md`](reproducibility/README.md) to run any of it.

## Install

```bash
pip install -e .
```

## Use

```python
from lntest import rank_genes_groups_ln

rank_genes_groups_ln(adata, groupby="leiden", layer="norm_counts")
adata.uns["rank_genes_groups"]["logfoldchanges"]
adata.uns["rank_genes_groups"]["lfc_se"]      # SE of the LFC; lfc ± 1.96·se is the interval
```

Takes **normalised, not log-transformed** counts.

## Tests

```bash
pytest
```
