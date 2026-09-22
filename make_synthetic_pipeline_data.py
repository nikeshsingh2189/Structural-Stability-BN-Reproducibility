"""Writes two SYNTHETIC files in the public inspection-file schema (columns JL.m, RD.m, SPD.m, EPD.m, SPO.deg, SIPRD.m,
SIPO.deg, L.mm, W.mm, MD.mm, SO.deg) drawn from the manuscript's fixed DAG.  For testing the scripts only: the published
numbers are NOT expected to reproduce on these files.   python scripts/make_synthetic_pipeline_data.py OUTDIR [N]"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tests"))
from bn_radius.pipeline_fit import ORDER, ANGULAR
from test_pipeline_fit import simulate_codes

out = sys.argv[1]; N = int(sys.argv[2]) if len(sys.argv) > 2 else 30000
os.makedirs(out, exist_ok=True)
states = [8 if v in ANGULAR else 4 for v in ORDER]
for year, seed in ((1, 11), (7, 17)):
    rng = np.random.default_rng(seed)
    codes = simulate_codes(N, rng, states); c = {v: codes[:, i] for i, v in enumerate(ORDER)}
    u = lambda: rng.random(N); JL = 12.0
    df = pd.DataFrame({"JL.m": JL, "RD.m": (c["UR"] + u()) / 4 * JL, "SPD.m": 1.0})
    df["EPD.m"] = df["SPD.m"] + (c["UE"] + u()) / 4 * JL - df["RD.m"]
    df["SPO.deg"] = c["SPO"] * 45 + 45 * u(); df["SIPRD.m"] = (c["US"] + u()) / 4 * JL; df["SIPO.deg"] = c["SIPO"] * 45 + 45 * u()
    df["L.mm"] = c["L"] + u(); df["W.mm"] = c["W"] + u(); df["MD.mm"] = c["MD"] + u()
    df["SO.deg"] = c["SO"] * 45 + 45 * u(); df.loc[rng.random(N) < 0.005, "SO.deg"] = np.nan
    df.to_csv(os.path.join(out, f"synthetic_year{year}.csv"), index=False)
print("wrote", out)
