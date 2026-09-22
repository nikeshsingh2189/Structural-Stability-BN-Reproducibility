"""Audit archived/generated non-timing outputs against manuscript/supplement values.

Timing values are deliberately excluded because wall-clock results are hardware dependent.
Raw-data pipeline checks are reported as SKIP unless the public Year-1/Year-7 files
have been supplied and the corresponding rerun outputs are present.
"""
from __future__ import annotations
import csv, os, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES = os.path.join(ROOT, "results")
checks = []

def add(label, ok, detail=""):
    checks.append((label, ok, detail))

def txt(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        return ""
    with open(p, encoding="utf-8") as f:
        return f.read()

# Exhaustive verification
s = txt("exhaustive_verification.txt")
for label, needle in [
    ("exhaustive total DAGs", "DAGs=33866"),
    ("exhaustive valid queries", "valid=1950439"),
    ("exhaustive mismatches", "mismatch=0 weighted_mismatch=0"),
    ("finite/infinite counts", "finite=945994 inf=1004445"),
    ("nonmonotone count", "nonmono=161416"),
    ("rho>=2 count", "rho>=2=53972"),
    ("six-node radius counts", "rho=1..4: 883486 51980 1522 16"),
]:
    add(label, needle in s)

# Larger random DAGs
s = txt("random_dag_check.txt")
add("random-DAG retained queries", "retained valid queries: 2000" in s)
add("random-DAG unit/weighted mismatches", "mismatches  unit-cost: 0   weighted: 0" in s)

# Benchmarks
s = txt("benchmarks_N10000.txt")
for label, needle in [
    ("Asia benchmark", "queries=665  rho=1:521 (78.35%)  rho=2:95 (14.29%)  rho=3:7 (1.05%)"),
    ("Asia nonmonotone", "nonmonotone:205 (30.83%)"),
    ("Insurance benchmark", "queries=10000  rho=1:8689 (86.89%)  rho=2:112 (1.12%)  rho=3:1 (0.01%)"),
    ("Insurance nonmonotone", "nonmonotone:1065 (10.65%)"),
    ("ALARM benchmark", "queries=10000  rho=1:2156 (21.56%)  rho=2:221 (2.21%)  rho=3:12 (0.12%)  rho>=4:3 (0.03%)  rho=inf:7608 (76.08%)"),
    ("ALARM nonmonotone", "nonmonotone:941 (9.41%)"),
]:
    add(label, needle in s)

# ALARM representative and all minimum sets
p = os.path.join(RES, "alarm_monitoring_cases.csv")
if os.path.exists(p):
    with open(p, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    by_pair = {(r["X"], r["Y"]): r for r in rows}
    add("ALARM displayed HRSAT representative",
        by_pair.get(("HRSAT", "VENTMACH"), {}).get("one_min_breaking_deletion") == "INTUBATION;SHUNT;VENTTUBE")
else:
    add("ALARM displayed HRSAT representative", False, "missing CSV")

p = os.path.join(RES, "alarm_all_minimum_breaking_sets.csv")
if os.path.exists(p):
    with open(p, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    by_pair = {(r["X"], r["Y"]): r for r in rows}
    add("ALARM HRSAT displayed set among all minima",
        "INTUBATION;SHUNT;VENTTUBE" in by_pair.get(("HRSAT", "VENTMACH"), {}).get("all_min_breaking_deletions", ""))
else:
    add("ALARM HRSAT displayed set among all minima", False, "missing CSV")

# ALARM protection
p = os.path.join(RES, "alarm_protection_design.csv")
if os.path.exists(p):
    with open(p, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    got = [round(float(r["exact_pct"]), 2) for r in rows[:5]]
    add("ALARM exact protection B=0..4", got == [87.96, 91.06, 95.10, 98.20, 100.00], str(got))
else:
    add("ALARM exact protection B=0..4", False, "missing CSV")

s = txt("alarm_protection_design.txt")
for p_loss, triple in [(0.1, "0 | 0 | 0"), (0.2, "2 | 7 | 7"), (0.3, "2 | 7 | 11"), (0.4, "3 | 7 | 13"), (0.5, "3 | 7 | 15")]:
    add(f"B95 p_loss={p_loss}", f"p_loss={p_loss}:" in s and f"B95 = {triple}" in s)

# Fixed pipeline DAG
s = txt("pipeline_radius_enumeration.txt")
add("pipeline total/radius counts", "valid relations: 1130   radius counts: {1: 1038, 2: 90, 3: 2}" in s)
add("pipeline SO baseline", "SO in C: 556   break after deleting SO: 190 (34.17%)   preserve: 366 (65.83%)" in s)

# Optional raw-data pipeline baseline fit
s = txt("pipeline_full_bn_run.txt")
if s:
    differs = "DIFFERS" in s
    required = [
        "edge agreement (Jaccard)", "valid relations", "relations with SO in C",
        "break after deleting SO", "preserved after deleting SO"
    ]
    add("raw pipeline baseline fit", (not differs) and all(x in s for x in required),
        "contains DIFFERS" if differs else "")
else:
    checks.append(("raw pipeline baseline fit", None, "not run; public Year-1/Year-7 files are external to this archive"))

# Optional raw-data sensitivity fit
s = txt("pipeline_sensitivity_run.txt")
if s:
    rule_ok = "prespecified DAG rule: identical Year-1/Year-7 fit if equal; otherwise their edge intersection" in s
    common_ok = all(x in s for x in [
        "common valid SO-conditioned relations: 556 (published 556)",
        "break under all three: 120 (published 120)",
        "under a subset: 134 (published 134)",
        "none: 302 (published 302)",
    ])
    add("raw pipeline sensitivity", rule_ok and common_ok and "DIFFERS" not in s,
        "prespecified intersection rule" if rule_ok else "prespecified-rule marker not found")
else:
    checks.append(("raw pipeline sensitivity", None, "not run; run with the public Year-1/Year-7 files"))

# Optional deterministic common-graph rerun
s = txt("pipeline_common_relation_sensitivity.txt")
if s:
    ok = all(x in s for x in [
        "common valid SO-conditioned relations :   556", "published 556   OK",
        "break under all three graphs        :   120", "published 120   OK",
        "break under a subset of the graphs  :   134", "published 134   OK",
        "survive under all three             :   302", "published 302   OK",
    ])
    add("common-556 fitted-graph analysis", ok)
else:
    checks.append(("common-556 fitted-graph analysis", None, "not run; generated after the public-data sensitivity fit"))

fails = [x for x in checks if x[1] is False]
for label, ok, detail in checks:
    tag = "PASS" if ok is True else ("SKIP" if ok is None else "FAIL")
    print(f"[{tag}] {label}" + (f" -- {detail}" if detail else ""))
print(f"\nSummary: {sum(x[1] is True for x in checks)} passed, {len(fails)} failed, {sum(x[1] is None for x in checks)} skipped")
if fails:
    raise SystemExit(1)
