#!/usr/bin/env python3
"""Fetch the pinned SMART-DS-derived P3R edge list and save source/target columns."""
from pathlib import Path
from urllib.request import urlopen
import io
import pandas as pd

COMMIT = "7f989a1f424e384dcd22d0ab10a76d4b5fce90ee"
URL = f"https://raw.githubusercontent.com/UCF-RISES/OutageMap/{COMMIT}/P3R/edgeList.csv"
OUT = Path(__file__).resolve().parent / "data" / "smartds_sfo_p3r_edge_pairs.csv"


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(URL, timeout=60) as response:
        raw = response.read()
    df = pd.read_csv(io.BytesIO(raw))
    if not {"source", "target"}.issubset(df.columns):
        raise RuntimeError("Pinned source file does not contain source/target columns")
    edges = df[["source", "target"]].copy()
    edges.to_csv(OUT, index=False)
    print(f"Wrote {len(edges)} directed edges to {OUT}")
    print(f"Unique nodes: {len(set(edges['source']).union(set(edges['target'])))}")

if __name__ == "__main__":
    main()
