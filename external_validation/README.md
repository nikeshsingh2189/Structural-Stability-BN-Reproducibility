# SMART-DS external structural validation package

This archive reproduces the external engineering-network validation reported for
the d-separation deletion-radius manuscript.

## Reported deterministic results

| Experiment | Relations | Exact agreement |
|---|---:|---:|
| Random valid relations, |C| = 1,...,6 | 1,000 | 1,000/1,000 |
| Stratified radii r = 1,2,3,4 (250 each) | 1,000 | 1,000/1,000 |
| Positive variable-specific deletion costs | 500 | 500/500 |
| **Total** | **2,500** | **2,500/2,500** |

The random set contains 927 radius-1, 68 radius-2, and 5 radius-3 relations.

## Graph

The processed SFO/P3R test graph contains 767 nodes and 766 directed edges. It
is a DAG and one weakly connected component. See `DATA_PROVENANCE.md` for the
public source and exact pinned Git identifiers.

## Files to archive with the paper

- `bn_radius/core.py`: exact core graph routines used in the manuscript.
- `run_external_validation_all.py`: recomputes all three experiments.
- `verify_reported_results.py`: independently recomputes all 2,500 stored cases
  and checks that the archived result CSVs are internally correct.
- `generate_random_queries.py`: regenerates the exact random 1,000-query set
  using seed `2026092201`.
- `data/smartds_sfo_p3r_edge_pairs.csv`: processed source-target topology.
- `data/*query_set*.csv`: fixed query inputs for exact reproducibility.
- `results/*.csv`: reported relation-level results and summaries.

## Reproduction

From this directory:

```bash
python -m pip install -r requirements.txt
python verify_reported_results.py
python generate_random_queries.py
python run_external_validation_all.py
```

Expected verification messages include:

```text
Verified 1000 random + 1000 stratified + 500 weighted cases: 2500/2500 exact agreement.
exact query-set match: True
All 2500 cases reproduced exactly.
```

## What is deterministic

The deletion radii/costs and agreement indicators are deterministic. Runtime
columns are hardware- and load-dependent and therefore will not reproduce bit
for bit on another machine.

## Interpretation

This package supports an **external structural validation** claim on a realistic
synthetic engineering-network topology. It should not be described as real
utility empirical validation, nor as the exact 10,090-node subnetwork used in a
previous IISE Transactions anomaly-detection study.
