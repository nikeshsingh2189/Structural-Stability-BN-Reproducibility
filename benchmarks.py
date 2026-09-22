"""Section 5.2 / Supplement C: deletion-radius distributions on Asia, Insurance and ALARM.
Asia: full enumeration (all pairs, all nonempty C).  Insurance/ALARM: sampling protocol of the paper --
draw an unordered pair uniformly, |C| uniformly in 1..8, C uniformly among subsets of that size, keep valid
statements until N are retained (duplicates are discarded; use --keep-duplicates to keep them).
Seeds 2026090801 (Insurance) and 2026090802 (ALARM) with Python's random.Random, drawing vertex indices in the
networks' original file order, reproduce the manuscript's retained queries exactly (Table S counts).
Usage: python scripts/benchmarks.py [--n 10000] [--keep-duplicates]"""
import argparse, itertools, math, os, random, sys
from collections import Counter, defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from bn_radius import load_network, d_separated, radius_moral_graph

ap = argparse.ArgumentParser(); ap.add_argument("--n", type=int, default=10000); ap.add_argument("--keep-duplicates", action="store_true")
args = ap.parse_args()


def nonmonotone(par, X, Y, C):
    C = list(C); k = len(C)
    conn = [not d_separated(par, X, Y, {C[i] for i in range(k) if not m >> i & 1}) for m in range(1 << k)]   # after LOSING set m
    return any(conn[r1] and not conn[r2] for r1 in range(1 << k) if conn[r1] for r2 in range(1 << k) if r2 & r1 == r1 and r2 != r1)


def wilson(x, n, z=1.96):
    p = x / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * (c - h), 100 * (c + h)


def summarize(label, rows):
    n = len(rows); cnt = Counter(r for r, _, _ in rows); nm = sum(m for _, _, m in rows)
    key = lambda r: "inf" if r == float("inf") else int(r)
    c = Counter(key(r) for r, _, _ in rows)
    ge4 = sum(v for k, v in c.items() if k != "inf" and k >= 4)
    print(f"{label}: queries={n}  rho=1:{c[1]} ({100*c[1]/n:.2f}%)  rho=2:{c[2]} ({100*c[2]/n:.2f}%)  rho=3:{c[3]} ({100*c[3]/n:.2f}%)  "
          f"rho>=4:{ge4} ({100*ge4/n:.2f}%)  rho=inf:{c['inf']} ({100*c['inf']/n:.2f}%)  nonmonotone:{nm} ({100*nm/n:.2f}%)")
    lo, hi = wilson(c[1], n); lo2, hi2 = wilson(c["inf"], n)
    print(f"      95% Wilson  rho=1: {lo:.2f}-{hi:.2f}   rho=inf: {lo2:.2f}-{hi2:.2f}")
    print(f"      share with rho<=3 (finite): {100*(c[1]+c[2]+c[3])/n:.2f}%")


# ---- Asia: exact
names, idx, par = load_network("asia", order="bif"); n = len(names); rows = []
for X, Y in itertools.combinations(range(n), 2):
    others = [v for v in range(n) if v not in (X, Y)]
    for k in range(1, len(others) + 1):
        for C in itertools.combinations(others, k):
            if d_separated(par, X, Y, C):
                rows.append((radius_moral_graph(par, X, Y, C), len(C), nonmonotone(par, X, Y, C)))
summarize("ASIA (exact enumeration)", rows)
print("   paper: 665 queries; rho=1:521 rho=2:95 rho=3:7 rho>=4:0 inf:42; nonmonotone 205 (30.83%)")

# ---- sampled networks
for name, seed in (("insurance", 2026090801), ("alarm", 2026090802)):
    names, idx, par = load_network(name, order="bif"); n = len(names); rng = random.Random(seed)   # original file order + Python RNG = the paper's sampler
    rows, seen = [], set()
    while len(rows) < args.n:
        X, Y = sorted(rng.sample(range(n), 2)); k = rng.randint(1, 8)
        C = tuple(sorted(rng.sample([v for v in range(n) if v not in (X, Y)], k)))
        if not d_separated(par, X, Y, C): continue
        if not args.keep_duplicates:
            if (X, Y, C) in seen: continue
            seen.add((X, Y, C))
        rows.append((radius_moral_graph(par, X, Y, C), k, nonmonotone(par, X, Y, C)))
    summarize(name.upper() + f" (sampled, N={args.n})", rows)
    strata = defaultdict(list)
    for r in rows: strata[r[1]].append(r)
    os.makedirs("results", exist_ok=True)
    with open(f"results/benchmark_strata_{name}.csv", "w", newline="") as fh:                # all |C| = 1..8 rows
        import csv; wr = csv.writer(fh); wr.writerow(["network", "|C|", "n", "rho=1_pct", "rho=2-3_pct", "rho>=4_finite_pct", "rho=inf_pct"])
        for k in range(1, 9):
            s = strata[k]; c = Counter(("inf" if r[0] == float("inf") else int(r[0])) for r in s); m = max(len(s), 1)
            ge4 = sum(v for kk, v in c.items() if kk != "inf" and kk >= 4)
            wr.writerow([name, k, len(s), f"{100*c[1]/m:.2f}", f"{100*(c[2]+c[3])/m:.2f}", f"{100*ge4/m:.2f}", f"{100*c['inf']/m:.2f}"])
    for k in (1, 4, 8):
        s = strata[k]; c = Counter(("inf" if r[0] == float("inf") else int(r[0])) for r in s)
        print(f"      |C|={k}: n={len(s)}  rho=1 {100*c[1]/len(s):.2f}%  rho=2-3 {100*(c[2]+c[3])/len(s):.2f}%  inf {100*c['inf']/len(s):.2f}%")
