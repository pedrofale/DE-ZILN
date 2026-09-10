# LN's $t$-test

[![CI](https://github.com/okviman/lntest/actions/workflows/ci.yml/badge.svg)](https://github.com/okviman/lntest/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/ln-ttest.svg)](https://pypi.org/project/ln-ttest/)

Differential expression testing on an asymptotically unbiased log-fold-change estimator.

## Installation

```bash
pip install ln-ttest
```

## Usage

```python
from lntest import rank_genes_groups_ln

rank_genes_groups_ln(adata, groupby="leiden", layer="norm_counts")
adata.uns["rank_genes_groups"]["logfoldchanges"]
adata.uns["rank_genes_groups"]["lfc_se"]      # SE of the LFC; lfc ± 1.96·se is the interval
```

Takes **normalised, not log-transformed** counts.
