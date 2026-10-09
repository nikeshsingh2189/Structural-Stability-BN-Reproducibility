#!/usr/bin/env python3
"""Recompute the SMART-DS SFO/P3R external structural validation.

The script uses the exact manuscript `bn_radius.core` routines bundled in this
archive and three pre-specified query sets. It recomputes all graph-theoretic
quantities from the public processed source-target edge list.

Timing columns are measured anew and are expected to vary by machine. Radius,
cost, and agreement columns are deterministic.
"""
from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

import networkx as nx
import pandas as pd

from bn_radius.core import (
    d_separated,
    deletion_radius,
    deletion_radius_bruteforce,
    radius_moral_graph,
)


def parse_C(value):
    text = "" if pd.isna(value) else str(value).strip()
    return tuple() if not text else tuple(int(x) for x in text.split(";"))


def parse_weights(value):
    out = {}
    for token in str(value).split(";"):
        k, v = token.split(":")
        out[int(k)] = float(v)
    return out


def same(a, b, tol=1e-10):
    if math.isinf(float(a)) and math.isinf(float(b)):
        return True
    return abs(float(a) - float(b)) <= tol


def load_graph(edge_file: Path):
    ed = pd.read_csv(edge_file)
    G = nx.DiGraph()
    G.add_edges_from(ed[["source", "target"]].itertuples(index=False, name=None))
    if not nx.is_directed_acyclic_graph(G):
        raise RuntimeError("Input network is not a DAG")
    nodes = sorted(G.nodes())
    par = {v: list(G.predecessors(v)) for v in nodes}
    return G, par


def graph_facts(G):
    return {
        "nodes": G.number_of_nodes(),
        "directed_edges": G.number_of_edges(),
        "is_DAG": nx.is_directed_acyclic_graph(G),
        "weak_components": nx.number_weakly_connected_components(G),
        "max_indegree": max(dict(G.in_degree()).values()),
        "max_outdegree": max(dict(G.out_degree()).values()),
        "sources": sum(d == 0 for _, d in G.in_degree()),
        "sinks": sum(d == 0 for _, d in G.out_degree()),
        "underlying_is_tree": nx.is_tree(G.to_undirected()),
    }


def run_unweighted(par, query_df, stratified=False):
    rows = []
    for qid, r in query_df.reset_index(drop=True).iterrows():
        x, y, C = int(r.X), int(r.Y), parse_C(r.C)
        if not d_separated(par, x, y, C):
            raise RuntimeError(f"Query {qid} is not initially d-separated")

        t0 = time.perf_counter()
        ra, Rstar = deletion_radius(par, x, y, C)
        t1 = time.perf_counter()
        rm = radius_moral_graph(par, x, y, C)
        t2 = time.perf_counter()
        rb, sets = deletion_radius_bruteforce(par, x, y, C)
        t3 = time.perf_counter()

        row = {
            "query_id": qid,
            "X": x,
            "Y": y,
            "C": ";".join(map(str, C)),
            "C_size": len(C),
            "radius_algorithm": ra,
            "radius_explicit_moral": rm,
            "radius_bruteforce": rb,
            "agree_all": same(ra, rm) and same(ra, rb),
            "algorithm_ms": (t1 - t0) * 1000,
            "explicit_moral_ms": (t2 - t1) * 1000,
            "bruteforce_ms": (t3 - t2) * 1000,
            "Rstar": ";".join(map(str, sorted(Rstar))) if Rstar is not None else "",
            "num_min_break_sets": len(sets),
        }
        if stratified:
            row["target_radius"] = int(r.target_radius)
            row["matches_target_radius"] = same(ra, int(r.target_radius))
        rows.append(row)
    return pd.DataFrame(rows)


def run_weighted(par, query_df):
    rows = []
    for _, r in query_df.reset_index(drop=True).iterrows():
        qid, x, y, C = int(r.query_id), int(r.X), int(r.Y), parse_C(r.C)
        weights = parse_weights(r.weights)
        if set(weights) != set(C):
            raise RuntimeError(f"Weight keys do not match C for query {qid}")
        if not d_separated(par, x, y, C):
            raise RuntimeError(f"Weighted query {qid} is not initially d-separated")

        t0 = time.perf_counter()
        ra, Rstar = deletion_radius(par, x, y, C, weights=weights)
        t1 = time.perf_counter()
        rm = radius_moral_graph(par, x, y, C, weights=weights)
        t2 = time.perf_counter()
        rb, sets = deletion_radius_bruteforce(par, x, y, C, weights=weights)
        t3 = time.perf_counter()
        rows.append({
            "query_id": qid,
            "X": x,
            "Y": y,
            "C": ";".join(map(str, C)),
            "C_size": len(C),
            "weights": r.weights,
            "cost_algorithm": ra,
            "cost_explicit_moral": rm,
            "cost_bruteforce": rb,
            "agree_all": same(ra, rm) and same(ra, rb),
            "algorithm_ms": (t1 - t0) * 1000,
            "explicit_moral_ms": (t2 - t1) * 1000,
            "bruteforce_ms": (t3 - t2) * 1000,
            "Rstar": ";".join(map(str, sorted(Rstar))) if Rstar is not None else "",
            "num_min_break_sets": len(sets),
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=Path(__file__).resolve().parent, type=Path)
    ap.add_argument("--outdir", default=None, type=Path)
    args = ap.parse_args()

    root = args.root.resolve()
    out = (args.outdir or (root / "reproduced_results")).resolve()
    out.mkdir(parents=True, exist_ok=True)

    G, par = load_graph(root / "data" / "smartds_sfo_p3r_edge_pairs.csv")
    pd.DataFrame([graph_facts(G)]).to_csv(out / "graph_facts_reproduced.csv", index=False)

    random_q = pd.read_csv(root / "data" / "random_query_set_1000.csv")
    strat_q = pd.read_csv(root / "data" / "stratified_query_set_1000.csv")
    weight_q = pd.read_csv(root / "data" / "weighted_query_set_500.csv")

    random_res = run_unweighted(par, random_q, stratified=False)
    strat_res = run_unweighted(par, strat_q, stratified=True)
    weight_res = run_weighted(par, weight_q)

    random_res.to_csv(out / "random_validation_reproduced.csv", index=False)
    strat_res.to_csv(out / "stratified_validation_reproduced.csv", index=False)
    weight_res.to_csv(out / "weighted_validation_reproduced.csv", index=False)

    summary = pd.DataFrame([
        ["Random valid relations", len(random_res), int(random_res.agree_all.sum())],
        ["Stratified deletion radii r=1,...,4", len(strat_res), int((strat_res.agree_all & strat_res.matches_target_radius).sum())],
        ["Variable-specific deletion costs", len(weight_res), int(weight_res.agree_all.sum())],
    ], columns=["experiment", "relations", "exact_agreements"])
    summary["agreement_rate"] = summary.exact_agreements / summary.relations
    summary.to_csv(out / "validation_summary_reproduced.csv", index=False)

    print(pd.DataFrame([graph_facts(G)]).to_string(index=False))
    print()
    print(summary.to_string(index=False))
    if int(summary.exact_agreements.sum()) != int(summary.relations.sum()):
        raise SystemExit("At least one external-validation case failed")
    print(f"\nAll {int(summary.relations.sum())} cases reproduced exactly.")


if __name__ == "__main__":
    main()
