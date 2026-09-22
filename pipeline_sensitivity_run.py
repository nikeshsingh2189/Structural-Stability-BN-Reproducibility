"""Section 5.4 / Supplement D: sensitivity of the SO-loss result to discretisation and
maximum indegree, plus the relation-level partition over parent caps 1, 2, and 3.

Primary reproduction command:

  python scripts/pipeline_sensitivity_run.py --year1 PATH --year7 PATH [--cut-side right|left]
      [--param-count possible|observed] [--filter-col COL --filter-val VAL] [--save-dags DIR]

Graph rule (prespecified, not selected from manuscript targets):
  * fit the ordered BIC DAG independently to Year 1 and Year 7;
  * if the two fitted edge sets are identical, use that common graph;
  * otherwise use their edge intersection as the stable descriptive DAG for that scenario.

This is the same stable-edge rule used by scripts/pipeline_full_bn_run.py.  Published values below
are used only after computation as validation checks; they never determine which graph is analysed.
"""
import argparse, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from bn_radius.pipeline_fit import (read_table, build_nodes, cut_points, discretise, fit_dag, jaccard,
                                    enumerate_relations, summarize, write_edge_csv)

ap = argparse.ArgumentParser()
ap.add_argument("--year1", required=True)
ap.add_argument("--year7", required=True)
ap.add_argument("--cut-side", default="right", choices=["right", "left"])
ap.add_argument("--param-count", default="possible", choices=["possible", "observed"])
ap.add_argument("--filter-col")
ap.add_argument("--filter-val")
ap.add_argument("--save-dags", metavar="DIR",
                help="write the prespecified quartile fitted graphs for caps 1-3 to DIR/pipeline_dags_cap{1,2,3}.csv")
a = ap.parse_args()

# Manuscript/Supplement values. These are validation targets only.
PUB = {
    (3, 2): (1.000, 1130, 556, 190, 34.17),
    (4, 1): (0.778, 2904, 1360, 288, 21.18),
    (4, 2): (1.000, 1130, 556, 190, 34.17),
    (4, 3): (0.882, 1002, 556, 254, 45.68),
    (5, 2): (1.000, 1130, 556, 190, 34.17),
}
PUB_PART = {"common": 556, "all": 120, "some": 134, "none": 302,
            "cap3_only": 64, "cap23": 70, "rates": (21.58, 34.17, 45.68)}


def stable_graph(e1, e7):
    """Prespecified scenario graph: common fit if identical, otherwise edge intersection."""
    return e1 if e1 == e7 else (e1 & e7)


raw = {y: read_table(p, a.filter_col, a.filter_val) for y, p in ((1, a.year1), (7, a.year7))}
nodes = {y: build_nodes(raw[y]) for y in raw}

fits, chosen, summaries = {}, {}, {}
for bins, cap in PUB:
    cuts = cut_points([nodes[1], nodes[7]], bins)
    e = {}
    for y in (1, 7):
        codes, states = discretise(nodes[y], cuts, bins, a.cut_side)
        e[y] = fit_dag(codes, states, cap, a.param_count)
    fits[(bins, cap)] = e
    chosen[(bins, cap)] = stable_graph(e[1], e[7])
    summaries[(bins, cap)] = summarize(enumerate_relations(chosen[(bins, cap)]))

print("prespecified DAG rule: identical Year-1/Year-7 fit if equal; otherwise their edge intersection")
print("bins cap  edge-agr.  chosen-arcs  valid  SO-panel  break  break%   | published row")
all_ok = True
for (bins, cap), p in PUB.items():
    e = fits[(bins, cap)]
    j = jaccard(e[1], e[7])
    s = summaries[(bins, cap)]
    pct = 100 * s["so_break"] / max(s["so_panel"], 1)
    ok = (round(j, 3) == p[0] and
          (s["valid"], s["so_panel"], s["so_break"]) == p[1:4] and
          abs(pct - p[4]) < 0.01)
    all_ok &= ok
    print(f"{bins:4d} {cap:3d}  {j:8.3f}  {len(chosen[(bins, cap)]):11d}"
          f"  {s['valid']:5d}  {s['so_panel']:8d}  {s['so_break']:5d}  {pct:6.2f}"
          f"   {'OK' if ok else 'DIFFERS'} | {p}")

if not all_ok:
    sys.exit("Prespecified sensitivity analysis differs from at least one archived manuscript value. "
             "No fitted DAGs were saved; verify preprocessing/data version before submission.")

if a.save_dags:
    os.makedirs(a.save_dags, exist_ok=True)
    for cap in (1, 2, 3):
        write_edge_csv(os.path.join(a.save_dags, f"pipeline_dags_cap{cap}.csv"), chosen[(4, cap)])
    print("saved prespecified fitted graphs to", a.save_dags)

# Figure 6: relations with SO in C that are valid under all three quartile parent-cap graphs.
rel = {cap: enumerate_relations(chosen[(4, cap)]) for cap in (1, 2, 3)}
panel = {cap: {k: v[1] for k, v in rel[cap].items() if v[1] is not None} for cap in rel}
common = set(panel[1]) & set(panel[2]) & set(panel[3])
pat = {}
for k in common:
    key = tuple(int(panel[cap][k]) for cap in (1, 2, 3))
    pat[key] = pat.get(key, 0) + 1
n = len(common)
allb = pat.get((1, 1, 1), 0)
none = pat.get((0, 0, 0), 0)
some = n - allb - none
rates = [100 * sum(panel[c][k] for k in common) / max(n, 1) for c in (1, 2, 3)]
part_ok = (n == PUB_PART["common"] and allb == PUB_PART["all"] and some == PUB_PART["some"] and
           none == PUB_PART["none"] and pat.get((0, 0, 1), 0) == PUB_PART["cap3_only"] and
           pat.get((0, 1, 1), 0) == PUB_PART["cap23"] and
           all(abs(g-r) < 0.01 for g, r in zip(rates, PUB_PART["rates"])))

print(f"\ncommon valid SO-conditioned relations: {n} (published {PUB_PART['common']})")
print(f"  break under all three: {allb} (published {PUB_PART['all']});"
      f"  under a subset: {some} (published {PUB_PART['some']});"
      f"  none: {none} (published {PUB_PART['none']})")
print(f"  patterns (cap1,cap2,cap3): {dict(sorted(pat.items()))}"
      f"   published partial patterns: cap3 only {PUB_PART['cap3_only']}, caps 2&3 {PUB_PART['cap23']}")
print(f"  break rate over the common relations, caps 1/2/3:"
      f" {rates[0]:.2f} / {rates[1]:.2f} / {rates[2]:.2f} %   published {PUB_PART['rates']}")
print("  note: the cap-1 rate over its full SO-panel is 288/1,360 = 21.18%; Figure 6 uses the common 556-relation set.")
if not part_ok:
    sys.exit("Prespecified common-relation sensitivity differs from the archived Figure 6 values.")
