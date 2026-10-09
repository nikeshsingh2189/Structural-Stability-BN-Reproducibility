"""Section 5.3: information protection on the five ALARM relations.
 * exact deletion-aware allocation (maximises mean preservation, all loss subsets enumerated)
 * radius-guided greedy heuristic (tie rule: most relations whose residual radius increases, then alphabetical)
 * degree-based rule (best value among admissible ties) and exact random-allocation mean
 * B95 budgets for p_loss in {0.1,...,0.5}; tie-resolution sensitivity (300 random tie resolutions)
 * protection stress design (Supplement C): radii after protecting VENTTUBE / VENTTUBE+CATECHOL
Writes results/alarm_protection_design.csv"""
import csv, itertools, math, os, random, sys
from functools import lru_cache
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from alarm_cases import load_cases
from bn_radius import d_separated, radius_moral_graph

names, idx, par, rel = load_cases()
items = sorted({v for _, _, C in rel for v in C}, key=lambda v: names[v])       # alphabetical order of variable names
N = len(items)
sep_tab = []
for X, Y, C in rel:
    k = len(C)
    sep_tab.append(np.array([d_separated(par, X, Y, set(C) - {C[i] for i in range(k) if m >> i & 1}) for m in range(1 << k)]))


def make_pi(p):
    @lru_cache(None)
    def pi_r(r, prot):
        C = rel[r][2]; k = len(C); free = [i for i in range(k) if C[i] not in prot]; tot = 0.0
        for sub in range(1 << len(free)):
            mask = 0; cnt = 0
            for j, i in enumerate(free):
                if sub >> j & 1: mask |= 1 << i; cnt += 1
            if sep_tab[r][mask]: tot += p ** cnt * (1 - p) ** (len(free) - cnt)
        return tot
    def Pi(S):
        S = set(S); return sum(pi_r(r, frozenset(S & set(rel[r][2]))) for r in range(len(rel))) / len(rel)
    return pi_r, Pi


def residual_radii(S):
    return [radius_moral_graph(par, X, Y, C, protected=S) for X, Y, C in rel]


def surrogate(S, p):                          # H(S) = sum_r p^{rho_r(S)} (Eq. 15), rho=inf contributes 0
    return sum(0.0 if r == float("inf") else p ** r for r in residual_radii(S))


def radius_guided(B, p, rng=None):
    S = []
    for _ in range(B):
        cand = [v for v in items if v not in S]
        sc = {v: surrogate(set(S) | {v}, p) for v in cand}; m = min(sc.values())
        ties = [v for v in cand if abs(sc[v] - m) < 1e-12]
        if rng is not None:
            S.append(rng.choice(ties)); continue                       # random tie resolution (sensitivity study)
        cur = residual_radii(set(S))
        gain = {v: sum(a != b for a, b in zip(residual_radii(set(S) | {v}), cur)) for v in ties}
        best = max(gain.values())
        S.append(min((v for v in ties if gain[v] == best), key=lambda v: names[v]))   # remaining ties: alphabetical by name
    return S


deg = {v: len(par[v]) + sum(v in par[w] for w in par) for v in items}
def degree_rule(B, Pi):
    if B == 0: return set()
    order = sorted(items, key=lambda v: (-deg[v], names[v])); cut = deg[order[B - 1]]
    must = [v for v in items if deg[v] > cut]; ties = [v for v in items if deg[v] == cut]
    return max((set(must) | set(t) for t in itertools.combinations(ties, B - len(must))), key=Pi)

def random_mean(B, pi_r):
    tot = 0.0
    for r in range(len(rel)):
        C = rel[r][2]; k = len(C); s = 0.0
        for t in range(k + 1):
            if not 0 <= B - t <= N - k: continue
            w = math.comb(k, t) * math.comb(N - k, B - t) / math.comb(N, B)
            s += w * np.mean([pi_r(r, frozenset(T)) for T in itertools.combinations(C, t)])
        tot += s
    return tot / len(rel)


os.makedirs("results", exist_ok=True)
rows = []
p = 0.30
pi_r, Pi = make_pi(p)
print(f"--- p_loss={p}: mean preservation (%)  [exact | radius-guided | degree | random mean]")
for B in range(0, 5):
    ex = max(itertools.combinations(items, B), key=lambda S: Pi(S)); rg = radius_guided(B, p)
    row = (B, [names[v] for v in ex], 100 * Pi(ex), [names[v] for v in rg], 100 * Pi(rg), 100 * Pi(degree_rule(B, Pi)), 100 * random_mean(B, pi_r))
    rows.append(row); print(f"B={B}: {row[2]:6.2f} | {row[4]:6.2f} | {row[5]:6.2f} | {row[6]:6.2f}   exact={row[1]}  radius-guided={row[3]}")
with open("results/alarm_protection_design.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["B", "exact_set", "exact_pct", "radius_guided_set", "radius_guided_pct", "degree_pct", "random_mean_pct"])
    for r in rows: w.writerow([r[0], ";".join(r[1]), f"{r[2]:.2f}", ";".join(r[3]), f"{r[4]:.2f}", f"{r[5]:.2f}", f"{r[6]:.2f}"])

print("\n--- B95 budgets (exact | degree | random mean), unprotected preservation")
for pl in (0.1, 0.2, 0.3, 0.4, 0.5):
    pr, P_ = make_pi(pl)
    b_ex = next(B for B in range(N + 1) if max(P_(S) for S in itertools.combinations(items, B)) >= 0.95 - 1e-12) if True else None
    b_dg = next(B for B in range(N + 1) if P_(degree_rule(B, P_)) >= 0.95 - 1e-12)
    b_rm = next(B for B in range(N + 1) if random_mean(B, pr) >= 0.95 - 1e-12)
    print(f"p_loss={pl}: unprotected {100*P_(set()):.2f}%  B95 = {b_ex} | {b_dg} | {b_rm}")

print("\n--- tie-resolution sensitivity of the radius-guided heuristic (300 random tie resolutions, p_loss=0.30)")
rng = random.Random(0)
for B in (2, 3, 4):
    vals = sorted(100 * Pi(set(radius_guided(B, p, rng))) for _ in range(300))
    print(f"B={B}: min {vals[0]:.2f}  median {vals[150]:.2f}  max {vals[-1]:.2f}  share at 100%: {sum(v >= 99.995 for v in vals)/300:.3f}")

print("\n--- protection stress design: exact residual radii (INT-HYPO, MINVOL-SAO2, CO-PRESS, HRSAT-VENT, HYPO-ERRLOW)")
for prot in ([], ["VENTTUBE"], ["VENTTUBE", "CATECHOL"]):
    S = {i for v in prot for i, n in enumerate(names) if n == v}
    print(f"protected={prot or 'none'}: {residual_radii(S)}")
