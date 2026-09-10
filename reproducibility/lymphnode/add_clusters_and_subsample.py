#!/usr/bin/env python3
"""
Post-processing of the Visium HD 2um-spot-to-segmented-cell assignment.

Step A (--out-annotated):
    Take the h5ad produced by `assign_visium_shapes.py --annotate-no-aggregate`
    (one row per 2um spot, obs['shape_id'] = segmented cell id) and attach the
    graph-based cluster label of the parent cell, taken from the GeoParquet
    written by `extract_cluster_from_geojson.py`.
    obs becomes ['cell_id', 'Cluster'] with a plain integer index.
    -> data/vishd-cluster1-cluster3-2um-with-clusters.h5ad

Step B (--out-subsampled):
    Keep only the spots whose parent cell has ID < --max-cell-id (default
    40000). Cell IDs in the 10x segmentation GeoJSON are handed out in raster
    order, so this is a spatially contiguous crop of the tissue (the top
    band of the section) and every kept cell keeps all of its 2um spots.
    -> data/vishd-cluster1-cluster3-2um-with-clusters-subsampled.h5ad
"""
import argparse

import anndata as ad
import geopandas as gpd
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-h5ad", required=True, help="h5ad with one row per 2um spot and obs['shape_id'].")
    p.add_argument("--shapes-parquet", required=True, help="GeoParquet with ['ID', 'Cluster'] (from extract_cluster_from_geojson.py).")
    p.add_argument("--out-annotated", required=True, help="Output h5ad with obs ['cell_id', 'Cluster'].")
    p.add_argument("--out-subsampled", required=True, help="Output h5ad restricted to cell_id < --max-cell-id.")
    p.add_argument("--max-cell-id", type=int, default=40000, help="Keep cells with ID strictly below this (default: 40000).")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    adata = ad.read_h5ad(args.in_h5ad)
    print(f"Loaded {adata.shape[0]} spots x {adata.shape[1]} genes from {args.in_h5ad}")

    shapes = gpd.read_parquet(args.shapes_parquet)
    cluster_of_cell = pd.Series(shapes["Cluster"].values, index=shapes["ID"].values)

    # Step A: spot -> parent cell id + cluster label, integer row index
    cell_id = pd.Index(adata.obs["shape_id"].to_numpy())
    obs = pd.DataFrame(
        {
            "cell_id": cell_id.to_numpy(),
            "Cluster": pd.Categorical(cluster_of_cell.reindex(cell_id).to_numpy()),
        }
    )
    if obs["Cluster"].isna().any():
        raise ValueError("Some spots have a shape_id that is absent from the shapes parquet.")
    adata.obs = obs
    print(adata.obs["Cluster"].value_counts().to_string())
    adata.write(args.out_annotated)
    print(f"Wrote annotated h5ad: {args.out_annotated}")

    # Step B: spatial crop by cell id
    keep = adata.obs["cell_id"].to_numpy() < args.max_cell_id
    sub = adata[keep].copy()
    # subsetting stringifies the integer row index; keep it integer as in the original
    sub.obs.index = sub.obs.index.astype("int64")
    print(
        f"Subsampled to cell_id < {args.max_cell_id}: {sub.shape[0]} spots, "
        f"{sub.obs['cell_id'].nunique()} cells"
    )
    print(sub.obs["Cluster"].value_counts().to_string())
    sub.write(args.out_subsampled)
    print(f"Wrote subsampled h5ad: {args.out_subsampled}")


if __name__ == "__main__":
    main()
