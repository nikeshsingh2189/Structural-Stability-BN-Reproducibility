"""Supplement F: exact random-order first-break counts (iid continuous loss times => all |C|! orders equally likely).
Writes results/alarm_random_loss_exact.csv"""
import csv, itertools, os, sys
from collections import Counter
from functools import lru_cache
sys.path.insert(0, os.path.dirname(__file__))
from alarm_cases import load_cases, CASES
from bn_radius import d_separated

names, idx, par, rel = load_cases()
os.makedirs("results", exist_ok=True)
with open("results/alarm_random_loss_exact.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["X", "Y", "orders", "first_break_counts", "E[L|L<inf]", "Pr(L=rho)_pct", "Pr(never)_pct"])
    for (x, y, _), (X, Y, C) in zip(CASES, rel):
        k = len(C)

        @lru_cache(None)
        def connected(lost_mask):
            return not d_separated(par, X, Y, {C[i] for i in range(k) if not lost_mask >> i & 1})

        cnt = Counter()
        for perm in itertools.permutations(range(k)):
            mask, L = 0, None
            for j, i in enumerate(perm, 1):
                mask |= 1 << i
                if connected(mask):
                    L = j
                    break
            cnt[L] += 1
        total = sum(cnt.values())
        fin = [(L, c) for L, c in cnt.items() if L is not None]
        EL = sum(L * c for L, c in fin) / sum(c for _, c in fin) if fin else float("nan")
        show = {("never" if L is None else L): c for L, c in sorted(cnt.items(), key=lambda t: (t[0] is None, t[0]))}
        p_rho = 100 * cnt[min(L for L, _ in fin)] / total if fin else float("nan")
        p_never = 100 * cnt[None] / total
        print(f"{x:>12}-{y:<13} orders={total:6d} {show}  E[L|L<inf]={EL:.2f}  Pr(L=rho)={p_rho:.2f}%  Pr(never)={p_never:.2f}%")
        w.writerow([x, y, total, show, f"{EL:.2f}", f"{p_rho:.2f}", f"{p_never:.2f}"])
