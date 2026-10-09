# Reproducibility archive

## Structural Stability of Conditional Independence in Bayesian Networks under Conditioning Information Loss

This archive contains the implementation, verification experiments, benchmark analyses, ALARM protection calculations, pipeline-analysis code, and SMART-DS external structural-validation materials associated with the manuscript.

## Quick start

```bash
python3 -m pip install -r requirements.txt
python3 -m pytest -q tests
bash run_all.sh
```

Without the external pipeline spreadsheets, `run_all.sh` reproduces the theorem checks, benchmark analyses, ALARM calculations, and fixed-DAG pipeline enumeration. The public pipeline fitting and graph-sensitivity analyses run when the Year-1 and Year-7 spreadsheets are supplied.

```bash
export YEAR1="/path/External Feature Year One.xlsx"
export YEAR7="/path/External Feature Year Seven.xlsx"
export YEAR3="/path/External Feature Year Three.xlsx"   # optional
export YEAR5="/path/External Feature Year Five.xlsx"    # optional
bash run_all.sh
```

The public pipeline data are available from Mendeley Data DOI `10.17632/c2h2jf5c54.1` under CC BY 4.0. The spreadsheets are not redistributed in this archive.

## SMART-DS external validation

The complete external structural-validation package is under `external_validation/`. It contains the processed SFO/P3R topology, the fixed random/stratified/weighted query sets, all 2,500 relation-level reported results, and scripts that independently recompute the stored results.

From the repository root:

```bash
python external_validation/verify_reported_results.py
python external_validation/generate_random_queries.py
python external_validation/run_external_validation_all.py
```

Expected deterministic checks include `2500/2500 exact agreement`, exact regeneration of the 1,000-query random set, and exact reproduction of all deletion radii/costs. Runtime columns are hardware-dependent. See `external_validation/README.md` and `external_validation/DATA_PROVENANCE.md` for provenance and interpretation.

## Main archive layout

| Path | Purpose |
|---|---|
| `bn_radius/core.py` | d-separation, exhaustive deletion radius, fixed-graph characterization, sparse auxiliary graph, Algorithm 1 |
| `bn_radius/networks.py` | benchmark-network edge-list utilities |
| `bn_radius/pipeline_fit.py` | pipeline node construction, discretization, BIC fitting, graph enumeration, common-relation comparison |
| `cpp/radius_tools.cpp` | exhaustive verification, chain scaling, and high-indegree stress test |
| `scripts/random_dag_check.py` | larger random-DAG verification |
| `scripts/benchmarks.py` | Asia, Insurance, and ALARM benchmark calculations |
| `scripts/alarm_monitoring_cases.py` | ALARM monitoring relations and minimum breaking deletions |
| `scripts/alarm_protection_design.py` | protection-allocation calculations |
| `scripts/alarm_random_loss_exact.py` | exact random-order loss analysis |
| `scripts/pipeline_radius_enumeration.py` | enumeration on the reported fixed pipeline DAG |
| `scripts/pipeline_full_bn_run.py` | baseline pipeline graph fitting from public data |
| `scripts/pipeline_sensitivity_run.py` | discretization/parent-cap sensitivity and fitted-DAG export |
| `scripts/pipeline_common_relation_sensitivity.py` | common-relation sensitivity analysis |
| `scripts/audit_reported_numbers.py` | checks archived non-timing outputs against reported values |
| `external_validation/` | SMART-DS external structural-validation materials |
| `reference/` | manuscript target values used only for validation/documentation |
| `results/` | archived machine-readable/text outputs |
| `tests/` | unit tests |

## Random seeds

- larger random-DAG verification: `20260908`
- Insurance benchmark sampling: `2026090801`
- ALARM benchmark sampling: `2026090802`
- SMART-DS random experiment seed documented in the external-validation materials: `2026092201`

## Software and timing

Direct dependencies are pinned in `requirements.txt`; `requirements-lock.txt` records the environment used for the archived Python results. Wall-clock timings depend on hardware and system load and are not expected to match exactly on another machine.

## Pipeline sensitivity

The baseline parent-cap-2 edge list is included as `data/pipeline_dags_cap2.csv`. Parent-cap-1 and parent-cap-3 graphs are regenerated from the public Year-1 and Year-7 spreadsheets using a prespecified stable-edge rule: use the common fit when the two yearly edge sets are identical; otherwise use their edge intersection. Manuscript values are used only as post-computation validation checks. Run:

```bash
python3 scripts/pipeline_sensitivity_run.py \
  --year1 "$YEAR1" --year7 "$YEAR7" --save-dags data
python3 scripts/pipeline_common_relation_sensitivity.py
```

The common-relation analysis uses the 556 relations valid under all three fitted graph specifications. Target values under `reference/` are used only for validation and are not fitted-graph inputs.
