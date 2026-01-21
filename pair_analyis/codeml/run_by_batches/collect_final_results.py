#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
from statsmodels.stats.multitest import multipletests


partial = sorted(Path("results_partial").glob("result_*.tsv"))
dfs = [pd.read_csv(f, sep="\t") for f in partial]

df = pd.concat(dfs, ignore_index=True)
df["p_fdr"] = multipletests(df["p_raw"], alpha=0.05, method="fdr_bh")[1]

df.to_csv("branch_site_final.tsv", sep="\t", index=False)
print("Final table written.")