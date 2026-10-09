"""Loading benchmark networks (bnlearn repository structures) and small helpers.

The structures (variables and arcs only, no parameters) of Asia, Insurance and ALARM are shipped in data/networks/ so that
nothing has to be downloaded and results do not depend on the installed pgmpy version.  pgmpy is only a fallback."""
from __future__ import annotations
import csv
import os
import warnings
from typing import Dict, List, Sequence, Tuple

_DATA = os.path.join(os.path.dirname(__file__), "..", "data", "networks")


def parents_from_edges(edges: Sequence[Tuple[str, str]], nodes: Sequence[str] = ()):
    """Build (names, index, parents) from a list of directed edges (u -> v); names sorted alphabetically."""
    names = sorted(set(nodes) | {u for u, _ in edges} | {v for _, v in edges})
    idx = {v: i for i, v in enumerate(names)}
    par: Dict[int, List[int]] = {i: [] for i in range(len(names))}
    for u, v in edges:
        par[idx[v]].append(idx[u])
    return names, idx, par


def _local_edges(name: str):
    path = os.path.join(_DATA, f"{name}.edges.csv")
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        rows = list(csv.reader(f))[1:]
    return [(r[0], r[1]) for r in rows]


def load_network(name: str, order: str = "alpha"):
    """Load 'asia', 'insurance' or 'alarm'. Returns (names, idx, parents).

    order="alpha": variables sorted alphabetically (default).
    order="bif":   original file order (bn_radius/node_orders.py); needed to reproduce the paper's sampled benchmark queries."""
    edges = _local_edges(name)
    if edges is None:                                            # fallback: pgmpy download (needs internet)
        warnings.filterwarnings("ignore")
        try:
            from pgmpy.utils import get_example_model
        except ImportError:
            from pgmpy.example_models import load_model as get_example_model
        model = get_example_model(name)
        edges = [(u, v) for u, v in model.edges()]
    if order == "bif":
        from .node_orders import BIF_ORDER
        names = list(BIF_ORDER[name])
        idx = {v: i for i, v in enumerate(names)}
        par: Dict[int, List[int]] = {i: [] for i in range(len(names))}
        for u, v in edges:
            par[idx[v]].append(idx[u])
        return names, idx, par
    from .node_orders import BIF_ORDER
    return parents_from_edges(edges, BIF_ORDER[name])
