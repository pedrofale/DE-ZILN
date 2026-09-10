"""LN's t-test: differential expression on an asymptotically unbiased log-fold-change estimator.

    from lntest import get_LN_lfcs, rank_genes_groups_ln

``rank_genes_groups_ln`` is imported eagerly, and can be, because the package
depends on numpy and scipy alone. It briefly did not: ``scanpy_wrapper`` imported
statsmodels for one call and pandas for nothing at all, which would have made a
plain re-export here break ``import lntest`` for anyone without those. That was
solved by deleting the dependencies rather than by deferring the import.
"""
from __future__ import annotations

from .ln_test import (
    TRIGAMMA_EXACT,
    TRIGAMMA_RECOMB25,
    get_LN_lfcs,
    get_LN_lfcs_sparse,
    trigamma_diff_int,
    trigamma_diff_recomb25,
)
from .scanpy_wrapper import rank_genes_groups_ln

try:  # pragma: no cover - absent only for a source tree with no installed metadata
    from importlib.metadata import PackageNotFoundError, version

    __version__ = version("lntest")
except (ImportError, PackageNotFoundError):  # pragma: no cover
    __version__ = "0.0.0.dev0"

__all__ = [
    "TRIGAMMA_EXACT",
    "TRIGAMMA_RECOMB25",
    "get_LN_lfcs",
    "get_LN_lfcs_sparse",
    "rank_genes_groups_ln",
    "trigamma_diff_int",
    "trigamma_diff_recomb25",
    "__version__",
]
