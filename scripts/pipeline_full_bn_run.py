"""Section 5.4 / Supplement D: baseline pipeline analysis on the public inspection files.

  python scripts/pipeline_full_bn_run.py --year1 PATH --year7 PATH [--year3 PATH --year5 PATH] [--check-epd]
        [--bins 4 --cap 2] [--cut-side right|left] [--param-count possible|observed]
        [--filter-col COL --filter-val VAL]

Steps: build the nine nodes; pooled Year-1/Year-7 quantile cut points (eight 45-degree sectors for angular nodes);
order-constrained BIC search fitted independently to Year 1 and Year 7; edge agreement; full enumeration of the deletion
radius on the fixed DAG; Year-3 seam-orientation loss.  Every published value is printed next to the value obtained here.

The data files are NOT redistributed here (public release: Khan et al., 2021; Yarveisy et al., 2021).  Options marked
OPTION in bn_radius/pipeline_fit.py are conventions the manuscript does not pin down; if a number differs, try the
alternatives (the sensitivity script lists them side by side)."""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from bn_radius.pipeline_fit import (ORDER, RAW, PAPER_EDGES, read_table, build_nodes, cut_points, discretise, fit_dag, jaccard,
                                    enumerate_relations, summarize)

ap = argparse.ArgumentParser()
ap.add_argument("--year1", required=True); ap.add_argument("--year7", required=True)
ap.add_argument("--year3"); ap.add_argument("--year5")
ap.add_argument("--bins", type=int, default=4); ap.add_argument("--cap", type=int, default=2)
ap.add_argument("--cut-side", default="right", choices=["right", "left"])
ap.add_argument("--param-count", default="possible", choices=["possible", "observed"])
ap.add_argument("--filter-col"); ap.add_argument("--filter-val"); ap.add_argument("--check-epd", action="store_true")
ap.add_argument("--sweep", action="store_true", help="try every (cut-side, param-count) convention and report which reproduces the Supplement D DAG, then exit")
a = ap.parse_args()

PUB = {"rows": {1: 735755, 3: 832681, 5: 759683, 7: 870301}, "so_missing": {1: 3604, 5: 4224, 7: 4801},
       "edge_agreement": 1.0, "n_edges": 14, "valid": 1130, "rho": {1: 1038, 2: 90, 3: 2}, "so_panel": 556, "so_break": 190, "so_preserve": 366,
       "epd_mad": {1: 9.3e-6, 7: 8.1e-12}}


def chk(label, got, exp, tol=0.0):
    ok = (abs(got - exp) <= tol) if isinstance(exp, (int, float)) else (got == exp)
    print(f"  {label:<44} got {got!s:<18} published {exp!s:<14} {'OK' if ok else 'DIFFERS'}")
    return ok


raw = {}
for y, p in ((1, a.year1), (7, a.year7), (3, a.year3), (5, a.year5)):
    if p: raw[y] = read_table(p, a.filter_col, a.filter_val)
print("== data")
for y, df in raw.items():
    chk(f"Year {y}: rows", len(df), PUB["rows"][y])
    miss = len(df) if "SO.deg" not in df else int(df["SO.deg"].isna().sum())
    if y in PUB["so_missing"]: chk(f"Year {y}: rows with SO.deg missing", miss, PUB["so_missing"][y])
    else: print(f"  Year {y}: SO.deg {'absent' if 'SO.deg' not in df else 'present'} (published: absent throughout Year 3)")

nodes = {y: build_nodes(raw[y]) for y in (1, 7)}
cuts = cut_points([nodes[1], nodes[7]], a.bins)
if a.sweep:
    print(f"== convention sweep (bins={a.bins}, cap={a.cap}): does the fit equal the 14-edge DAG of Supplement D in both years?")
    for side in ("right", "left"):
        for pc in ("possible", "observed"):
            e = {}
            for y in (1, 7):
                codes, states = discretise(nodes[y], cuts, a.bins, side); e[y] = fit_dag(codes, states, a.cap, pc)
            print(f"  cut-side={side:<5} param-count={pc:<8}  Year1==paper: {e[1] == PAPER_EDGES!s:<5}  Year7==paper: {e[7] == PAPER_EDGES!s:<5}"
                  f"  |Year1 xor paper|={len(e[1] ^ PAPER_EDGES)}  |Year7 xor paper|={len(e[7] ^ PAPER_EDGES)}  Jaccard(Y1,Y7)={jaccard(e[1], e[7]):.3f}")
    sys.exit(0)
fit = {}
for y in (1, 7):
    codes, states = discretise(nodes[y], cuts, a.bins, a.cut_side)
    fit[y] = fit_dag(codes, states, a.cap, a.param_count)
print(f"\n== fitted DAGs (bins={a.bins}, cap={a.cap}, cut-side={a.cut_side}, param-count={a.param_count})")
for y in (1, 7): print(f"  Year {y}: {len(fit[y])} edges: {sorted(fit[y])}")
jac = jaccard(fit[1], fit[7])
chk("edge agreement (Jaccard)", round(jac, 3), PUB["edge_agreement"], 1e-9)
chk("number of edges", len(fit[1]), PUB["n_edges"])
print("  identical to the DAG printed in Supplement D:", fit[1] == PAPER_EDGES and fit[7] == PAPER_EDGES)
if fit[1] != PAPER_EDGES: print("  edges differing from Supplement D (Year 1):", sorted(fit[1] ^ PAPER_EDGES))

edges = fit[1] if fit[1] == fit[7] else fit[1] & fit[7]
print("\n== full enumeration on the " + ("common" if fit[1] == fit[7] else "intersection of Year-1/Year-7") + " DAG")
s = summarize(enumerate_relations(edges))
chk("valid relations", s["valid"], PUB["valid"])
for r in (1, 2, 3): chk(f"radius {r}", s["rho"].get(r, 0), PUB["rho"][r])
chk("radius infinity", s["rho"].get("inf", 0), 0)
chk("relations with SO in C", s["so_panel"], PUB["so_panel"])
chk("break after deleting SO", s["so_break"], PUB["so_break"])
chk("preserved after deleting SO", s["so_preserve"], PUB["so_preserve"])
print(f"  SO-loss break rate: {100 * s['so_break'] / max(s['so_panel'], 1):.2f}% (published 34.17%)")

if a.check_epd:
    print("\n== EPD.m - SPD.m versus L.mm/1000 (why the Year-5 EPD.m omission is not used)")
    for y in (1, 7):
        d = (raw[y]["EPD.m"] - raw[y]["SPD.m"] - raw[y]["L.mm"] / 1000.0).abs()
        chk(f"Year {y}: mean absolute discrepancy (m)", float(d.mean()), PUB["epd_mad"][y], PUB["epd_mad"][y] * 0.15)
        print(f"     Year {y}: max discrepancy {d.max():.6f} m (published Year 1: 0.001)")
