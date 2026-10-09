"""Figure 6 / Section 5.4: SO-loss sensitivity across the parent-cap-1, -2 and -3 fitted graphs on the common valid relations.

Deterministic: it needs only the three fitted edge lists (no raw data).  By default it reads
    data/pipeline_dags_cap1.csv, data/pipeline_dags_cap2.csv, data/pipeline_dags_cap3.csv      (header: parent,child)
cap 2 is the DAG printed in Supplement D.  The cap-1 and cap-3 lists are the fitted graphs used for Figure 6; write them with
    python scripts/pipeline_sensitivity_run.py --year1 ... --year7 ... --save-dags data
Outputs: results/pipeline_common_relation_sensitivity.csv (one row per common relation) and a printed/published-value comparison."""
import argparse, csv, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from bn_radius.pipeline_fit import read_edge_csv, common_relation_partition

ap = argparse.ArgumentParser()
ap.add_argument("--dags", nargs=3, metavar=("CAP1", "CAP2", "CAP3"),
                default=[f"data/pipeline_dags_cap{c}.csv" for c in (1, 2, 3)])
ap.add_argument("--out", default="results/pipeline_common_relation_sensitivity.csv")
a = ap.parse_args()

missing = [p for p in a.dags if not os.path.exists(p)]
if missing:
    sys.exit("missing fitted edge list(s): " + ", ".join(missing) +
             "\nCreate them with pipeline_sensitivity_run.py --save-dags data (needs the public inspection files).")
dags = [read_edge_csv(p) for p in a.dags]
rows, s = common_relation_partition(dags)

os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
with open(a.out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["X", "Y", "C", "breaks_cap1", "breaks_cap2", "breaks_cap3"])
    for X, Y, C, pat in rows: w.writerow([X, Y, ";".join(C), *pat])

PUB = {"common": 556, "all": 120, "some": 134, "none": 302, "rates": (21.58, 34.17, 45.68),
       "full": [(2904, 1360, 288), None, (1002, 556, 254)]}
ok = lambda g, e: "OK" if g == e else "DIFFERS"
print("graphs:", ", ".join(f"cap {c}: {len(e)} arcs" for c, e in zip((1, 2, 3), dags)))
print(f"common valid SO-conditioned relations : {s['common']:5d}   published {PUB['common']}   {ok(s['common'], PUB['common'])}")
print(f"  break under all three graphs        : {s['all']:5d}   published {PUB['all']}   {ok(s['all'], PUB['all'])}")
print(f"  break under a subset of the graphs  : {s['some']:5d}   published {PUB['some']}   {ok(s['some'], PUB['some'])}")
print(f"  survive under all three             : {s['none']:5d}   published {PUB['none']}   {ok(s['none'], PUB['none'])}")
print(f"  patterns (cap1, cap2, cap3) -> count: {dict(sorted(s['patterns'].items()))}   (published partial patterns: cap 3 only 64, caps 2&3 70)")
print("  break rate over the common relations: " + " / ".join(f"{r:.2f}" for r in s["rates"]) + f" %   published {PUB['rates']}")
for c, fp, pp in zip((1, 2, 3), s["full_panels"], PUB["full"]):
    line = f"  cap {c} full graph: valid {fp['valid']}, SO-panel {fp['so_panel']}, breaks {fp['so_break']} ({100*fp['so_break']/max(fp['so_panel'],1):.2f}%)"
    print(line + (f"   published {pp}  {ok((fp['valid'], fp['so_panel'], fp['so_break']), pp)}" if pp else ""))
print("wrote", a.out)
