#!/usr/bin/env python3
"""Measure what substituting ``lntest``'s estimator for ``utils.py``'s does to each arm.

``utils.py`` and ``lntest`` are not two copies of one estimator. They differ in
the trigamma term of the standard error -- ``utils.trigamma(x) = 1/x`` against
``lntest.trigamma_diff_int(a, n) = sum_{j=a..n-1} 1/j^2`` -- and every published
RECOMB number came from the approximate one. The log-fold change is untouched by
that difference; the standard error, the t statistic and every p-value are not.

This harness answers the only question that gates the substitution: **does any
arm's headline number move?** It does not fix anything, and it must not be made
to pass. See decision-retire-utils-py, and handoff-math-questions questions 1-2.

Two gates
---------
**Hard.** ``max|dLFC| < 1e-5`` on every measured arm. The LFC involves no
trigamma, so the only difference expected there is ``utils``'s float32
intermediates. Anything larger means something other than the trigamma term
changed, and the run is a failure, not a finding.

**Soft.** Significance-call disagreements at alpha = 0.05 are *expected* to be
nonzero. They are reported, never failed on. A clean zero here would be evidence
this harness is wrong, not that the substitution is safe.

**If a headline number moves, stop.** That is a manuscript change and it goes to
Oskar. Do not soften, round, or average a difference away to get past this.

Usage
-----
    python check_utils_vs_lntest.py --arm all --json report.json

Requires ``lntest`` importable (``pip install -e .`` at the repo root, or the
``de-ziln-reproducibility`` conda env). The comparison is deliberately against
the *installed* package, not ``pkg/``, because the published package is what the
refactor puts behind every arm.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import pathlib
import sys
from dataclasses import dataclass, asdict

import numpy as np

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent

# The hard gate on the log-fold change. Not a tolerance to be relaxed: see the
# module docstring for why exceeding it is a failure rather than a finding.
HARD_DLFC_TOL = 1e-5

ALPHA = 0.05

# The five functions utils_frozen.py vendors, and which this harness re-verifies
# against utils.py for as long as utils.py exists.
FROZEN_FUNCS = ["trigamma", "get_DELN_lfcs", "get_LN_lfcs",
                "compute_p_vals", "get_t_statistic"]

# Paths decision-repo-restructure deletes. Used to assert -- rather than grep and
# hope -- that the utils API surface with no lntest counterpart has no caller
# that will survive the refactor.
DELETION_LIST = (
    "notebooks/scrnaseq_reanalyses/",
    "notebooks/test/de_test_with_scores.py",
    "notebooks/test/de_celltype_vs_rest.py",
    "notebooks/test/de_cluster_comparison_gsea.py",
    "reproducibility/celltype/",
    "generate_lfc_data.py",
    "scanpy_confidence_interval_test.py",
)


# --------------------------------------------------------------------------
# Arm registry
# --------------------------------------------------------------------------

MEASURABLE = "measurable"
# decision-retire-utils-py requires these two strings verbatim in the record.
NOT_APPLICABLE = "not applicable"
COVERED_BY_ARGUMENT = "covered by argument, not measured"


@dataclass
class Arm:
    letter: str
    name: str
    script: str
    status: str
    # Which utils entry point the arm actually calls. Arm C is the one that does
    # not call get_LN_lfcs, which is why it is an estimator substitution rather
    # than an import rename.
    utils_entry: str | None
    correction: str | None
    headline: str | None
    # Decimal places the manuscript prints this arm's headline metric to. The
    # "did it move" test is applied at this precision, and the raw delta is
    # reported regardless so the judgement is reviewable.
    headline_decimals: int | None
    note: str = ""


# NOTE: the letters A-H are used throughout the decision records but are never
# defined in one place. This mapping is inferred from the records that do name an
# arm (C = CITE-seq, F = lymph node, G = clustering, H = celltype) plus the
# handoff's "the two kidney experiments" and "the synthetic-NB experiments".
# It is the harness's own statement of scope and should be confirmed.
ARMS = [
    Arm("A", "synthetic NB null (FPR vs dispersion)", "variance_vs_metric.py",
        MEASURABLE, "get_LN_lfcs", "bonferroni", "fpr", 2,
        note="Pure null: mu1 == mu2, so every gene is a true negative."),
    Arm("B", "synthetic NB with planted DEGs", "large_scale_NB_DE_test.py",
        MEASURABLE, "get_LN_lfcs", "bonferroni", "tpr,f1", 2,
        note="60 runs (3 dispersions x 20 reps). Headline is the per-dispersion mean."),
    Arm("C", "CITE-seq memory CD4", "large_scale_CITE_seq_exp.py",
        MEASURABLE, "get_DELN_lfcs", "fdr_bh", "tpr,f1", 3,
        note="Estimator substitution, not an import rename: calls get_DELN_lfcs, "
             "which has no lntest counterpart. Its eps=1e-9 term and missing "
             "all-zero guard must be shown to be inert behind the arm's own "
             "min_cells_per_group=3 filter, not assumed to be."),
    Arm("D", "Visium HD kidney, UMI subsampling",
        "notebooks/test/vishd_test_parallel.py", COVERED_BY_ARGUMENT,
        "get_LN_lfcs", None, None, None,
        note="Aliases get_LN_lfcs as get_DELN_lfcs -- the same function arms A "
             "and B measure. Gated on a multi-gigabyte 10x download that has "
             "not happened; no Visium HD data is present in the tree."),
    Arm("E", "Visium HD kidney, spot subsampling",
        "notebooks/test/vishd_test_shape_split_shared.py", COVERED_BY_ARGUMENT,
        "get_LN_lfcs", None, None, None,
        note="Same alias as arm D, same missing download."),
    Arm("F", "lymph node subsampling",
        "reproducibility/run_subsampling_analysis_parallel.py", NOT_APPLICABLE,
        None, None, None, None,
        note="Reaches the estimator through de_utils.run_ln_de -> "
             "rank_genes_groups_ln -> get_LN_lfcs_sparse and never imports "
             "utils. This is a proof of no-effect, not an untested gap."),
    Arm("G", "PBMC3k clustering",
        "reproducibility/clustering/run_clustering_de.py", NOT_APPLICABLE,
        None, None, None, None,
        note="Same proof as arm F. Run anyway for continuity with "
             "decision-reconcile-dense-sparse, and because it is the only arm "
             "with committed reference checksums."),
]


# --------------------------------------------------------------------------
# Self-checks: the harness must not be trusted more than it has earned
# --------------------------------------------------------------------------

def _extract_funcs(path: pathlib.Path, names):
    src = path.read_text()
    tree, lines = ast.parse(src), src.splitlines(keepends=True)
    return {n.name: "".join(lines[n.lineno - 1:n.end_lineno]).rstrip()
            for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names}


def verify_frozen_copy():
    """Confirm utils_frozen.py is still a verbatim copy of utils.py.

    After stage 3d deletes utils.py this can no longer be checked, and saying so
    is the point: a referee running this on the post-refactor tree should be told
    that the frozen copy is now the only witness, not left to assume it was
    verified.
    """
    utils = REPO_ROOT / "utils.py"
    frozen = pathlib.Path(__file__).with_name("utils_frozen.py")
    if not utils.exists():
        return ("unverifiable", "utils.py is gone (expected after stage 3d); "
                                "utils_frozen.py is now the only copy of the "
                                "RECOMB-era estimator and cannot be re-checked.")
    a, b = _extract_funcs(utils, FROZEN_FUNCS), _extract_funcs(frozen, FROZEN_FUNCS)
    drifted = [n for n in FROZEN_FUNCS if a.get(n) != b.get(n)]
    if drifted:
        raise SystemExit(
            f"FATAL: utils_frozen.py has drifted from utils.py in {drifted}. "
            "The frozen copy must be verbatim or every number below is measured "
            "against the wrong estimator."
        )
    digest = hashlib.sha256(utils.read_bytes()).hexdigest()
    return ("verified", f"all {len(FROZEN_FUNCS)} functions verbatim against "
                        f"utils.py sha256 {digest[:16]}")


def require_lntest():
    try:
        import lntest
    except ImportError as exc:
        raise SystemExit(
            f"FATAL: cannot import lntest ({exc}).\n"
            "This harness compares against the published package on purpose.\n"
            "Install it with `pip install -e .` at the repo root, or activate "
            "the de-ziln-reproducibility conda env."
        )
    return lntest


def check_api_mismatches(lntest):
    """Assert the three API mismatches rather than trusting a grep.

    decision-retire-utils-py lists these as the rewiring's real work. Each is
    asserted here so that a change to either module fails this harness rather
    than surfacing as a silent behaviour change during stage 3b.
    """
    findings = []

    # (1) return_log_abs_statistic has no counterpart in lntest, and its only
    #     callers must all be on the deletion list.
    callers = []
    for path in sorted(REPO_ROOT.rglob("*.py")):
        rel = path.relative_to(REPO_ROOT).as_posix()
        if rel in ("utils.py", "reproducibility/utils_frozen.py") or ".git/" in rel:
            continue
        if "return_log_abs_statistic=True" in path.read_text():
            callers.append(rel)
    survivors = [c for c in callers if not any(c.startswith(d) for d in DELETION_LIST)]
    findings.append({
        "mismatch": "return_log_abs_statistic has no lntest counterpart",
        "callers": callers,
        "callers_surviving_the_refactor": survivors,
        "ok": not survivors,
        "detail": "Nothing to translate if every caller is deleted."
                  if not survivors else
                  "A caller survives the refactor and has no lntest equivalent.",
    })

    # (2) return_standard_error and return_statistic are mutually exclusive in
    #     lntest: it returns on the first and never looks at the second, where
    #     utils returned a five-tuple. A guard that raises is owed here; adding it
    #     is a later stage, so this only records the current behaviour.
    rng = np.random.default_rng(0)
    tiny_Y = rng.integers(0, 20, size=(30, 8)).astype(float)
    tiny_X = rng.integers(0, 20, size=(30, 8)).astype(float)
    both = lntest.get_LN_lfcs(tiny_Y, tiny_X, return_standard_error=True,
                              return_statistic=True)
    findings.append({
        "mismatch": "return_standard_error and return_statistic are mutually exclusive",
        "returned_tuple_length": len(both),
        "ok": len(both) == 3,
        "detail": "lntest silently drops the statistic when both flags are set. "
                  "A guard that raises is owed (decision-retire-utils-py, "
                  "Consequences 2); it is not added here because stage 1b moves "
                  "and changes nothing.",
    })

    # (3) get_DELN_lfcs has no counterpart at all -- this is what makes arm C a
    #     substitution rather than a rename.
    findings.append({
        "mismatch": "get_DELN_lfcs has no lntest counterpart",
        "ok": not hasattr(lntest, "get_DELN_lfcs"),
        "detail": "Arm C must be measured, not assumed equivalent.",
    })
    return findings


# --------------------------------------------------------------------------
# The comparison itself
# --------------------------------------------------------------------------

def compare_estimators(lfc_u, p_u, lfc_l, p_l, correction, alpha=ALPHA):
    """Compare one paired run of the two estimators on identical inputs.

    ``_u`` is utils_frozen, ``_l`` is lntest. Both must have been given exactly
    the same X and Y; producing that is each arm adapter's job, and is the part
    of this harness most able to be silently wrong.
    """
    import statsmodels.stats.multitest as smm

    lfc_u, p_u = np.asarray(lfc_u, float), np.asarray(p_u, float)
    lfc_l, p_l = np.asarray(lfc_l, float), np.asarray(p_l, float)

    dlfc = np.abs(lfc_u - lfc_l)
    finite_dlfc = dlfc[np.isfinite(dlfc)]

    # p can be exactly 0, so log10 is only defined on part of the vector. Report
    # how many genes were excluded rather than flooring them into the statistic.
    both_pos = (p_u > 0) & (p_l > 0) & np.isfinite(p_u) & np.isfinite(p_l)
    dlog10p = np.abs(np.log10(p_u[both_pos]) - np.log10(p_l[both_pos]))

    # Each arm applies its own correction; significance is judged after it.
    def adjust(p):
        ok = np.isfinite(p)
        out = np.full_like(p, np.nan)
        if ok.any():
            out[ok] = smm.multipletests(p[ok], alpha=alpha, method=correction)[1]
        return out

    adj_u, adj_l = adjust(p_u), adjust(p_l)
    sig_u, sig_l = adj_u < alpha, adj_l < alpha

    return {
        "n_genes": int(lfc_u.size),
        "max_abs_dlfc": float(finite_dlfc.max()) if finite_dlfc.size else None,
        "median_abs_dlfc": float(np.median(finite_dlfc)) if finite_dlfc.size else None,
        "max_abs_dlog10p": float(dlog10p.max()) if dlog10p.size else None,
        "n_excluded_from_dlog10p": int((~both_pos).sum()),
        # Split by direction: which estimator called the gene significant.
        "disagree_utils_only": int((sig_u & ~sig_l).sum()),
        "disagree_lntest_only": int((sig_l & ~sig_u).sum()),
        "disagree_total": int((sig_u != sig_l).sum()),
        "nan_lfc_utils": int(np.isnan(lfc_u).sum()),
        "nan_lfc_lntest": int(np.isnan(lfc_l).sum()),
        "nan_p_utils": int(np.isnan(p_u).sum()),
        "nan_p_lntest": int(np.isnan(p_l).sum()),
    }


def headline_moved(value_u, value_l, decimals):
    """Did the published number change at the precision the manuscript prints?

    Reported alongside the raw delta, never instead of it.
    """
    return round(float(value_u), decimals) != round(float(value_l), decimals)


# --------------------------------------------------------------------------
# Arm adapters -- stage 1b-ii
# --------------------------------------------------------------------------

def run_arm(arm: Arm, lntest, reps=None):
    """Reproduce one arm's inputs, run both estimators, recompute its headline.

    Not implemented yet. Each adapter has to transcribe its script's input
    generation exactly -- same seed, same draw order, same filtering -- because
    the scripts are top-level code that calls plt.show() and writes to hardcoded
    paths, so they cannot be imported. That transcription is the single most
    likely place for this harness to be quietly wrong, which is why it is a
    separate reviewed commit rather than bundled in with the machinery above.
    """
    raise NotImplementedError(
        f"arm {arm.letter} adapter is stage 1b-ii; the machinery above is 1b-i"
    )


# --------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--arm", default="all",
                    help="arm letter (A-G) or 'all' (default)")
    ap.add_argument("--reps", type=int, default=None,
                    help="override replicate count; the headline metric is only "
                         "comparable to the manuscript at the full count, so a "
                         "reduced run is a smoke test and is labelled as one")
    ap.add_argument("--json", type=pathlib.Path, default=None,
                    help="write the full report here")
    args = ap.parse_args(argv)

    report = {"frozen_copy": None, "api_mismatches": None, "arms": {}}

    state, detail = verify_frozen_copy()
    report["frozen_copy"] = {"state": state, "detail": detail}
    print(f"[frozen copy] {state}: {detail}\n")

    lntest = require_lntest()
    report["api_mismatches"] = check_api_mismatches(lntest)
    for f in report["api_mismatches"]:
        print(f"[api] {'ok  ' if f['ok'] else 'FAIL'} {f['mismatch']}")
        print(f"      {f['detail']}")
    print()

    selected = [a for a in ARMS if args.arm == "all" or a.letter == args.arm.upper()]
    if not selected:
        raise SystemExit(f"no such arm: {args.arm}")

    unimplemented = []
    for arm in selected:
        entry = asdict(arm)
        if arm.status != MEASURABLE:
            # These two statuses are proofs and arguments, not gaps, and the
            # decision record requires them to be recorded in these words.
            print(f"[arm {arm.letter}] {arm.status}: {arm.name}")
            print(f"           {arm.note}")
            entry["result"] = None
        else:
            try:
                entry["result"] = run_arm(arm, lntest, reps=args.reps)
            except NotImplementedError as exc:
                print(f"[arm {arm.letter}] UNIMPLEMENTED: {exc}")
                entry["result"] = None
                unimplemented.append(arm.letter)
        report["arms"][arm.letter] = entry
    print()

    if args.json:
        args.json.write_text(json.dumps(report, indent=2))
        print(f"report written to {args.json}")

    if unimplemented:
        # Refuse to look like a pass. A harness that exits 0 with nothing
        # measured is worse than one that fails.
        print(f"\nINCOMPLETE: arms {', '.join(unimplemented)} not measured. "
              "No gate has been evaluated; this is not a pass.")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
