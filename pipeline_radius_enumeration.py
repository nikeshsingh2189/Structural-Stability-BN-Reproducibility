"""Section 5.4 / Supplement D: full deletion-radius enumeration on the fixed pipeline DAG (14 edges,
identical in the Year-1 and Year-7 fits).  This script does NOT fit the DAG from the inspection data; the
fitting scripts (order-constrained BIC search, discretisation) must be run on the public files."""
import itertools, os, sys
from collections import Counter
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from bn_radius import parents_from_edges, d_separated, deletion_radius_bruteforce, deletion_radius

EDGES = [("SO", "UR"), ("SO", "SPO"), ("SO", "L"), ("SPO", "L"), ("SPO", "W"), ("L", "W"), ("UR", "UE"),
         ("L", "UE"), ("UR", "US"), ("UE", "US"), ("SPO", "SIPO"), ("W", "SIPO"), ("L", "MD"), ("W", "MD")]
NODES = ["SO", "UR", "SPO", "L", "W", "UE", "US", "SIPO", "MD"]
names, idx, par = parents_from_edges(EDGES, NODES); n = len(names); so = idx["SO"]
cnt, total, so_q, so_break = Counter(), 0, 0, 0
for X, Y in itertools.combinations(range(n), 2):
    others = [v for v in range(n) if v not in (X, Y)]
    for k in range(1, len(others) + 1):
        for C in itertools.combinations(others, k):
            if d_separated(par, X, Y, C):
                total += 1; rho, R = deletion_radius(par, X, Y, C); cnt[rho] += 1
                if so in C:
                    so_q += 1; so_break += not d_separated(par, X, Y, set(C) - {so})
print(f"valid relations: {total}   radius counts: { {int(k) if k != float('inf') else 'inf': v for k, v in sorted(cnt.items())} }")
print(f"SO in C: {so_q}   break after deleting SO: {so_break} ({100*so_break/so_q:.2f}%)   preserve: {so_q-so_break} ({100*(so_q-so_break)/so_q:.2f}%)")
for x, y, C in [("UR", "SPO", ["SO"]), ("UR", "W", ["SO", "SPO", "L"]), ("US", "SIPO", ["SO", "UR", "SPO", "L", "W", "UE"])]:
    rho, sets = deletion_radius_bruteforce(par, idx[x], idx[y], [idx[c] for c in C])
    print(f"{x} _||_ {y} | {C}: rho={int(rho)}  all minimum breaking deletions: {[sorted(names[v] for v in s) for s in sorted(sets, key=sorted)]}")
