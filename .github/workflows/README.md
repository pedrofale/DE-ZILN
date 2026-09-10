# Releasing `ln-ttest`

**Two names, on purpose.** The PyPI distribution is `ln-ttest`; the import is
`lntest`, and so is the GitHub repository. PyPI refused `lntest` as "too similar
to an existing project" — `intest` exists, and PyPI treats `l`, `i` and `1` as
confusable to block typosquatting. So `pip install ln-ttest` gives you
`from lntest import ...`, the way scikit-learn imports as sklearn.

Where each one belongs: the *distribution* name appears in `pyproject.toml`'s
`[project] name`, in `importlib.metadata.version(...)` inside `__init__.py`, in
the publish environment URL, and in the publisher table below. Everywhere else —
the wheel-contents guard, the floor job's import — is the *import* name and must
stay `lntest`.

`release.yml` publishes to PyPI when a version tag is pushed:

```bash
git tag v0.1.0 && git push upstream v0.1.0
```

It builds, checks, tests, and only then publishes. `workflow_dispatch` runs
everything except the publish step, so a release can be rehearsed without
creating a tag.

## One-time setup, before the first release

Publishing uses **PyPI Trusted Publishing (OIDC)**, so no API token is stored in
this repository. Nothing works until a *pending publisher* is registered.

1. On PyPI: **Your projects → Publishing → Add a new pending publisher**

   | Field | Value |
   |---|---|
   | PyPI project name | `ln-ttest` |
   | Owner | `okviman` |
   | Repository name | `lntest` |
   | Workflow name | `release.yml` |
   | Environment name | `pypi` |

2. In this repository: **Settings → Environments → New environment → `pypi`**.
   Adding required reviewers there turns every release into an approval step,
   which is worth doing while the version numbers are still cheap.

## What the workflow checks, and why each check exists

- **`twine check --strict`** — catches metadata PyPI would reject *after* the
  tag exists, when it is too late to fix quietly.
- **tag matches `pyproject.toml`** — the tag is the release's only claim about
  its version and pyproject is the artifact's only claim. A mismatch puts the
  wrong number on PyPI permanently.
- **wheel contains exactly three modules** — setuptools reuses `build/lib` when
  it is present. On 2026-09-10 a wheel built after renaming `ln_test.py` and
  `scanpy_wrapper.py` shipped the old names *alongside* the new ones, leaving
  `lntest.ln_test` importable and masking the rename. The workflow deletes
  `build/` first and then asserts the contents anyway.
- **the `floor` job** — installs the package alone, without the `test` extra, on
  the exact version `requires-python` advertises, and uses it.

  It exists because the floor was wrong once. pyproject claimed `>=3.9` while a
  signature used `bool | None` — PEP 604, which evaluates at definition time —
  so the package was unimportable on the version it advertised. Found by reading
  on 2026-09-11, not by CI, because there was no CI. The floor is now `>=3.10`,
  which is what the syntax needs and no restriction for this audience: scanpy
  requires `>=3.10` and anndata `>=3.11`.

## What this does not do

There is **no CI on push or pull request** — the test matrix runs only on tags
and manual dispatch. Between releases nothing checks a change. A `ci.yml` that
runs the same `test` job on `push` and `pull_request` would close that, and is
the obvious next addition.
