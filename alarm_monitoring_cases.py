"""Table 3 / Supplement C: radius, one and ALL minimum breaking deletions, minimum moral path.
Writes results/alarm_monitoring_cases.csv and results/alarm_all_minimum_breaking_sets.csv"""
import csv, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from alarm_cases import load_cases, CASES
from bn_radius import deletion_radius, deletion_radius_bruteforce, d_separated

names, idx, par, rel = load_cases()
os.makedirs("results", exist_ok=True)
with open("results/alarm_monitoring_cases.csv", "w", newline="") as f1, \
     open("results/alarm_all_minimum_breaking_sets.csv", "w", newline="") as f2:
    w1, w2 = csv.writer(f1), csv.writer(f2)
    w1.writerow(["X", "Y", "C", "rho", "one_min_breaking_deletion"])
    w2.writerow(["X", "Y", "rho", "all_min_breaking_deletions"])
    for (x, y, Cn), (X, Y, C) in zip(CASES, rel):
        assert d_separated(par, X, Y, C)
        rho, R = deletion_radius(par, X, Y, C)                       # Algorithm 1 radius
        rho_bf, allmin = deletion_radius_bruteforce(par, X, Y, C)    # direct enumeration
        assert rho == rho_bf
        # The theorem can return several equally valid minimum breaking deletions.
        # For manuscript display, use a deterministic alphabetically first representative
        # so that the CSV matches the representative set printed in the paper.
        ordered = sorted(allmin, key=lambda s: sorted(names[v] for v in s))
        R_display = ordered[0] if ordered else None
        one = "" if R_display is None else ";".join(sorted(names[v] for v in R_display))
        sets = " | ".join(";".join(sorted(names[v] for v in s)) for s in ordered)
        w1.writerow([x, y, ";".join(Cn), rho, one]); w2.writerow([x, y, rho, sets])
        print(f"{x:>12} {y:<13} rho={rho}  one R*={one or '-'}   all min: {sets or 'none'}")
