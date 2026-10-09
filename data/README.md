# Pipeline fitted-graph files

`pipeline_dags_cap2.csv` is the baseline 14-edge DAG printed in Supplementary Material D.

The cap-1 and cap-3 fitted edge lists are generated from the public Year-1 and Year-7 spreadsheets rather than copied from manuscript target counts. Run:

```bash
python3 scripts/pipeline_sensitivity_run.py --year1 PATH_TO_YEAR1 --year7 PATH_TO_YEAR7 --save-dags data
```

For each scenario the script uses the common Year-1/Year-7 fit when the fits are identical and otherwise their edge intersection. This rule is applied before comparison with manuscript values. If the resulting analysis passes the archived validation checks, the script writes `pipeline_dags_cap1.csv`, `pipeline_dags_cap2.csv`, and `pipeline_dags_cap3.csv`. Then run `scripts/pipeline_common_relation_sensitivity.py` to regenerate Figure 6's common-relation classification.
