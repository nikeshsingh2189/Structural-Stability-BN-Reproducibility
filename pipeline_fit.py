"""Pipeline case study (Section 5.4, Supplement D): node construction, discretisation and order-constrained
BIC structure search, plus enumeration of the deletion radius on the fitted DAG.

Conventions that can affect raw-data reconstruction are exposed as explicit options (marked OPTION).
The primary reproduction scripts use fixed defaults and compare results with manuscript values only after computation."""
from __future__ import annotations
import itertools
import math
from typing import Dict, FrozenSet, List, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from .core import d_separated, deletion_radius
from .networks import parents_from_edges

ORDER = ["SO", "UR", "SPO", "L", "W", "UE", "US", "SIPO", "MD"]          # prespecified variable order
ANGULAR = {"SO", "SPO", "SIPO"}                                            # eight 45-degree sectors
RAW = ["JL.m", "RD.m", "SPD.m", "EPD.m", "SPO.deg", "SIPRD.m", "SIPO.deg", "L.mm", "W.mm", "MD.mm", "SO.deg"]
Edge = Tuple[str, str]                                                     # (parent, child)

# fixed descriptive DAG reported in the manuscript (Year-1 and Year-7 baseline fits are identical)
PAPER_EDGES: Set[Edge] = {("SO", "UR"), ("SO", "SPO"), ("SO", "L"), ("SPO", "L"), ("SPO", "W"), ("L", "W"), ("UR", "UE"),
                          ("L", "UE"), ("UR", "US"), ("UE", "US"), ("SPO", "SIPO"), ("W", "SIPO"), ("L", "MD"), ("W", "MD")}


# ------------------------------------------------------------------ data
def read_table(path: str, filter_col: str = None, filter_val: str = None) -> pd.DataFrame:
    df = pd.read_excel(path) if str(path).lower().endswith((".xlsx", ".xls")) else pd.read_csv(path)
    if filter_col:                                                         # OPTION: keep only e.g. external anomalies
        df = df[df[filter_col].astype(str) == str(filter_val)]
    return df.reset_index(drop=True)


def build_nodes(raw: pd.DataFrame) -> pd.DataFrame:
    """Nine descriptive nodes (Supplement D, Table S-pipenodes). Missing SO.deg (Year 3) gives an all-NaN SO."""
    g = lambda c: raw[c].astype(float) if c in raw else pd.Series(np.nan, index=raw.index)
    out = pd.DataFrame(index=raw.index)
    out["SO"] = g("SO.deg")
    out["UR"] = g("RD.m") / g("JL.m")
    out["SPO"] = g("SPO.deg")
    out["L"] = g("L.mm")
    out["W"] = g("W.mm")
    out["UE"] = (g("EPD.m") - g("SPD.m") + g("RD.m")) / g("JL.m")
    out["US"] = g("SIPRD.m") / g("JL.m")
    out["SIPO"] = g("SIPO.deg")
    out["MD"] = g("MD.mm")
    return out[ORDER].replace([np.inf, -np.inf], np.nan)


def cut_points(nodes_list: Sequence[pd.DataFrame], bins: int) -> Dict[str, np.ndarray]:
    """Pooled (Year-1 + Year-7) quantile cut points for the continuous nodes."""
    pooled = pd.concat(nodes_list, ignore_index=True)
    qs = [k / bins for k in range(1, bins)]
    return {v: np.quantile(pooled[v].dropna().to_numpy(), qs) for v in ORDER if v not in ANGULAR}


def discretise(nodes: pd.DataFrame, cuts: Dict[str, np.ndarray], bins: int, side: str = "right") -> Tuple[np.ndarray, List[int]]:
    """Integer codes (-1 = missing) and the number of declared states per node.
    OPTION side: 'right' -> code = #cuts <= x ; 'left' -> code = #cuts < x."""
    codes = np.full((len(nodes), len(ORDER)), -1, dtype=np.int64)
    states = []
    for j, v in enumerate(ORDER):
        x = nodes[v].to_numpy(dtype=float); ok = ~np.isnan(x)
        if v in ANGULAR:
            codes[ok, j] = np.floor((x[ok] % 360.0) / 45.0).astype(np.int64) % 8; states.append(8)
        else:
            codes[ok, j] = np.searchsorted(cuts[v], x[ok], side=side); states.append(bins)
    return codes, states


# ------------------------------------------------------------------ structure search
def local_bic(codes: np.ndarray, states: Sequence[int], child: int, parents: Sequence[int],
              param_count: str = "possible") -> float:
    """Discrete BIC of the family (child | parents) on the complete cases of that family.
    OPTION param_count: 'possible' -> (r_child - 1) * prod(r_parent) with declared state counts;
                        'observed' -> counts only the observed states / parent configurations."""
    cols = [child] + list(parents)
    sub = codes[:, cols]; ok = (sub >= 0).all(axis=1); sub = sub[ok]; N = len(sub)
    if N == 0:
        return -math.inf
    r = states[child]; q = 1; cfg = np.zeros(N, dtype=np.int64)
    for p in parents:
        cfg = cfg * states[p] + sub[:, 1 + list(parents).index(p)]; q *= states[p]
    joint = cfg * r + sub[:, 0]
    counts = np.bincount(joint, minlength=q * r).reshape(q, r).astype(float)
    nj = counts.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        ll = np.nansum(np.where(counts > 0, counts * np.log(counts / nj), 0.0))
    if param_count == "possible":
        k = (r - 1) * q
    else:
        k = (int((counts.sum(axis=0) > 0).sum()) - 1) * int((nj[:, 0] > 0).sum())
    return ll - 0.5 * k * math.log(N)


def fit_dag(codes: np.ndarray, states: Sequence[int], cap: int, param_count: str = "possible") -> Set[Edge]:
    """Order-constrained search: for each child, all parent subsets among earlier variables with size <= cap are scored
    and the BIC-maximising family is kept (ties: fewer parents first, then order of the variables)."""
    edges: Set[Edge] = set()
    for c in range(len(ORDER)):
        best, best_set = -math.inf, ()
        for k in range(0, min(cap, c) + 1):
            for P in itertools.combinations(range(c), k):
                s = local_bic(codes, states, c, P, param_count)
                if s > best + 1e-12:
                    best, best_set = s, P
        edges |= {(ORDER[p], ORDER[c]) for p in best_set}
    return edges


def jaccard(a: Set[Edge], b: Set[Edge]) -> float:
    return len(a & b) / len(a | b) if (a | b) else 1.0


# ------------------------------------------------------------------ deletion-radius enumeration on a DAG
def enumerate_relations(edges: Set[Edge]):
    """All valid relations X _||_ Y | C (X<Y in ORDER, C nonempty). Returns dict key -> (radius, breaks_after_deleting_SO)
    where the second entry is None if SO not in C. key = (X, Y, frozenset(C)) with node names."""
    names, idx, par = parents_from_edges(sorted(edges), ORDER)
    inv = {i: n for n, i in idx.items()}
    out = {}
    for X, Y in itertools.combinations(sorted(range(len(names)), key=lambda i: ORDER.index(inv[i])), 2):
        others = [v for v in range(len(names)) if v not in (X, Y)]
        for k in range(1, len(others) + 1):
            for C in itertools.combinations(others, k):
                if d_separated(par, X, Y, C):
                    rho, _ = deletion_radius(par, X, Y, C, check_valid=False)
                    so = None
                    if idx["SO"] in C:
                        so = not d_separated(par, X, Y, set(C) - {idx["SO"]})
                    out[(inv[X], inv[Y], frozenset(inv[c] for c in C))] = (rho, so)
    return out


def summarize(rel) -> dict:
    from collections import Counter
    cnt = Counter(int(r) if r != float("inf") else "inf" for r, _ in rel.values())
    so = [b for _, b in rel.values() if b is not None]
    return {"valid": len(rel), "rho": dict(cnt), "so_panel": len(so), "so_break": sum(so), "so_preserve": len(so) - sum(so)}


# ------------------------------------------------------------------ common-relation comparison (Figure 6)
def read_edge_csv(path: str) -> Set[Edge]:
    """Edge list with header 'parent,child' (one arc per row)."""
    df = pd.read_csv(path)
    return {(str(p), str(c)) for p, c in zip(df["parent"], df["child"])}


def write_edge_csv(path: str, edges: Set[Edge]) -> None:
    pd.DataFrame(sorted(edges), columns=["parent", "child"]).to_csv(path, index=False)


def common_relation_partition(dags: Sequence[Set[Edge]]):
    """Relations X _||_ Y | C with SO in C that are d-separated under EVERY given DAG (the 'common valid relations'),
    classified by the pattern of SO-loss breaks across the DAGs.

    Returns (rows, summary): rows = [(X, Y, C, pattern_tuple)], summary = dict with the counts."""
    panels = []
    full = []
    for e in dags:
        rel = enumerate_relations(e)
        full.append(summarize(rel))
        panels.append({k: v[1] for k, v in rel.items() if v[1] is not None})      # relation -> breaks after deleting SO
    common = set(panels[0])
    for p in panels[1:]:
        common &= set(p)
    rows, pat = [], {}
    for k in sorted(common, key=lambda k: (k[0], k[1], sorted(k[2]))):
        pattern = tuple(int(p[k]) for p in panels)
        rows.append((k[0], k[1], sorted(k[2]), pattern))
        pat[pattern] = pat.get(pattern, 0) + 1
    m = len(dags)
    n = len(common)
    all_break = pat.get((1,) * m, 0)
    none_break = pat.get((0,) * m, 0)
    summary = {"common": n, "all": all_break, "none": none_break, "some": n - all_break - none_break, "patterns": pat,
               "rates": [100.0 * sum(p[k] for k in common) / n if n else float("nan") for p in panels], "full_panels": full}
    return rows, summary
