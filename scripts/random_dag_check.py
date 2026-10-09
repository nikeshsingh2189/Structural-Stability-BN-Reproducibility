"""Supplement B: larger-graph check. 50 random DAGs with 8-20 vertices (several densities), 2,000 retained
valid queries with |C| <= 8, unit-cost and integer-weighted (0..15) comparisons against direct deletion.
Seed 20260908 as in the supplement (Python's RNG differs from the authors', so the draws differ)."""
import os, random, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from bn_radius import d_separated, deletion_radius, deletion_radius_bruteforce

rng = random.Random(20260908)
retained = mism_unit = mism_w = bad_R = 0
for g in range(50):
    n = rng.randint(8, 20); dens = rng.choice([0.08, 0.12, 0.2, 0.3])
    par = {v: [u for u in range(v) if rng.random() < dens] for v in range(n)}
    perm = list(range(n)); rng.shuffle(perm)                       # hide the topological order
    par = {perm[v]: [perm[u] for u in par[v]] for v in par}
    kept, tries = 0, 0
    while kept < 40 and tries < 20000:
        tries += 1
        X, Y = rng.sample(range(n), 2); rest = [v for v in range(n) if v not in (X, Y)]
        C = rng.sample(rest, rng.randint(1, min(8, len(rest))))
        if not d_separated(par, X, Y, C): continue
        kept += 1; retained += 1
        w = {c: rng.randint(0, 15) for c in C}
        for weights, tag in ((None, "unit"), (w, "weighted")):
            brute, _ = deletion_radius_bruteforce(par, X, Y, C, weights)
            rho, R = deletion_radius(par, X, Y, C, weights)
            if rho != brute:
                mism_unit += tag == "unit"; mism_w += tag == "weighted"
            elif R is not None and d_separated(par, X, Y, set(C) - set(R)):
                bad_R += 1
print(f"retained valid queries: {retained}")
print(f"mismatches  unit-cost: {mism_unit}   weighted: {mism_w}   returned deletion not breaking: {bad_R}")
