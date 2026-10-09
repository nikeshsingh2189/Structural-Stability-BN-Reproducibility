"""Mechanical tests of the pipeline module on synthetic data drawn from the manuscript's fixed DAG."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from bn_radius.pipeline_fit import (ORDER, ANGULAR, PAPER_EDGES, build_nodes, cut_points, discretise, fit_dag, jaccard,
                                    enumerate_relations, summarize)


def simulate_codes(n, rng, states, edges=PAPER_EDGES, conc=0.25):
    """Sample node codes from a random discrete BN with the given edges (variables in ORDER)."""
    parents = {v: [p for p, c in sorted(edges) if c == v] for v in ORDER}
    data = {}
    for v in ORDER:
        r = states[ORDER.index(v)]; ps = parents[v]
        q = int(np.prod([states[ORDER.index(p)] for p in ps])) if ps else 1
        cpt = rng.dirichlet(np.full(r, conc), size=q)
        cfg = np.zeros(n, dtype=np.int64)
        for p in ps:
            cfg = cfg * states[ORDER.index(p)] + data[p]
        u = rng.random(n)
        data[v] = (u[:, None] > np.cumsum(cpt[cfg], axis=1)).sum(axis=1).clip(0, r - 1)
    return np.column_stack([data[v] for v in ORDER])


def test_paper_dag_enumeration():
    s = summarize(enumerate_relations(PAPER_EDGES))
    assert s["valid"] == 1130 and s["rho"] == {1: 1038, 2: 90, 3: 2}
    assert (s["so_panel"], s["so_break"], s["so_preserve"]) == (556, 190, 366)


def test_search_recovers_known_dag():
    rng = np.random.default_rng(0)
    states = [8 if v in ANGULAR else 4 for v in ORDER]
    codes = simulate_codes(60000, rng, states)
    # a little missingness, handled by complete-case scoring
    miss = rng.random(codes.shape) < 0.01; codes = np.where(miss, -1, codes)
    edges = fit_dag(codes, states, cap=2)
    assert jaccard(edges, PAPER_EDGES) >= 0.9, (edges ^ PAPER_EDGES)


def test_cap_one_gives_sparser_graph():
    rng = np.random.default_rng(1)
    states = [8 if v in ANGULAR else 4 for v in ORDER]
    codes = simulate_codes(40000, rng, states)
    e1, e2 = fit_dag(codes, states, cap=1), fit_dag(codes, states, cap=2)
    assert len(e1) <= 8 and len(e1) < len(e2)                      # at most one parent per child (8 non-root nodes)


def test_node_mapping_and_discretisation():
    raw = pd.DataFrame({"JL.m": [12.0] * 4, "RD.m": [3.0, 6.0, 9.0, 12.0], "SPD.m": [1.0] * 4, "EPD.m": [2.0, 3.0, 4.0, 5.0],
                        "SPO.deg": [10, 50, 100, 359.0], "SIPRD.m": [6.0] * 4, "SIPO.deg": [0, 44.9, 45.0, 90.0],
                        "L.mm": [1.0, 2.0, 3.0, 4.0], "W.mm": [1.0, 2.0, 3.0, 4.0], "MD.mm": [1.0, 2.0, 3.0, 4.0]})
    nodes = build_nodes(raw)
    assert list(nodes.columns) == ORDER and nodes["SO"].isna().all()               # SO.deg absent => all missing
    assert np.allclose(nodes["UR"], [0.25, 0.5, 0.75, 1.0]) and np.allclose(nodes["UE"], (raw["EPD.m"] - 1 + raw["RD.m"]) / 12)
    codes, states = discretise(nodes, cut_points([nodes], 4), 4)
    assert list(codes[:, ORDER.index("SIPO")]) == [0, 0, 1, 2]                    # 45-degree sectors
    assert list(codes[:, ORDER.index("SPO")]) == [0, 1, 2, 7]
    assert (codes[:, 0] == -1).all() and states[0] == 8
