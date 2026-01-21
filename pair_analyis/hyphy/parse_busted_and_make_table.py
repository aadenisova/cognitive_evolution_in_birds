import json
import glob
import pandas as pd
import os
from statsmodels.stats.multitest import multipletests

# files = glob.glob("absrel_results_correct_inno_10/*.json")
files = glob.glob("busted_results_correct_inno/*.json")

rows = []

for f in files:
    gene = os.path.basename(f).replace(".json", "")

    if os.path.getsize(f) == 0:
        continue

    with open(f) as fh:
        data = json.load(fh)


    LRT = data["test results"]["LRT"]
    p_val = data["test results"]["p-value"]

    fits = data['fits']["Unconstrained model"]['Rate Distributions']

    rates_background = fits['Background']
    rates_test = fits['Test']
    rates_syn = fits["Synonymous site-to-site rates"]

    
    error_sink_syn = sum(
        v["proportion"]
        for v in rates_syn.values()
        if v["rate"] > 5   # или >10, зависит от твоей калибровки
    )
    
    possel_test_prop = 0
    possel_test_maxomega = 0
    for v in rates_test.values():
        if v["omega"] > 1:
            possel_test_prop += v["proportion"]
            possel_test_maxomega = max(possel_test_maxomega, v["omega"])

    omega0_test_prop = sum(
        v["proportion"]
        for v in rates_test.values()
        if v["omega"] == 0
    )

    omega_0_proportion = 0
    for cat, values in rates_syn.items():
        omega = values.get('rate', 0)
        if omega == 0:  
            omega_0_proportion += values.get('proportion', 0)

    error_sink_proportion_test = 0
    for cat, values in rates_test.items():
        omega = values.get('omega', 0)
        if omega >= 100:  # error-sink
            error_sink_proportion_test += values.get('proportion', 0)

    rows.append({
            "gene": gene,
            "p-value": p_val,
            "LRT": LRT,
            "error_sink_back": error_sink_proportion_test,
            "omega_0_proportion": omega_0_proportion,
            "error_sink_test": error_sink_proportion_test,
        })

df = pd.DataFrame(rows)

df = df[
    #(df["LRT"] > 5) 
    #& (df["p-value"] < 0.05)
    #& (df["possel_test_prop"] > 0.01)
    #& (df["max_omega_test"] > 5)
    (df["error_sink"] == 0)
    #& (df["omega0_test_prop"] < 0.95)
]
 
# ---- FDR внутри каждой ВЕТВИ (вида) по всем генам ----
df["fdr"] = multipletests(df["p-value"], method="fdr_bh")[1]
# df["fdr"] = df["p-value"].apply((lambda p: multipletests(p, method="fdr_bh")[1]))

df.sort_values(by = ["fdr"]).sort_values(
    by=["error_sink_back", "error_sink_test"], 
    ascending=True
).to_csv("busted_e_full_table.csv", index=False)
# df_filtered.to_csv("busted_e_genes_all_branches_selected.csv", index=False)

# print("Готово.")


# import json
# import glob
# import pandas as pd
# import os
# from statsmodels.stats.multitest import multipletests


# # RESULT_DIR = "/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/absrel_results"
# files = glob.glob(f"absrel_results_correct/*.json", )

# print(files)

# rows = []

# for f in files:
#     gene = f.split("/")[-1].replace(".json", "")

#     if os.path.getsize(f) == 0:
#         continue
#     with open(f) as fh:
#         data = json.load(open(f))

#     # список тестируемых ветвей → test
#     tested = data["tested"]["0"]
#     test_branches = [b for b, t in tested.items() if t == "test"]

#     branch_data = data["branch attributes"]["0"]

#     for branch in test_branches:
#         pval = branch_data[branch]["Uncorrected P-value"]
#         rows.append({
#             "gene": gene,
#             "branch": branch,
#             "raw_p": pval
#         })
       

# # создаём таблицу
# df = pd.DataFrame(rows)

# # FDR
# df["fdr"] = (
#     df.groupby("branch")["raw_p"]
#       .transform(lambda p: multipletests(p, method="fdr_bh")[1])
# )

# # df["fdr"] = multipletests(df["raw_p"], method="fdr_bh")[1]

# # Фильтр: во всех ветвях гена FDR < 0.05
# significant_genes = (
#     df.groupby("gene")
#       .apply(lambda x: (x["fdr"] < 0.05).all())
# )

# sig_genes = significant_genes[significant_genes].index.tolist()

# summary = (
#     df.assign(sig = df["fdr"] < 0.05)
#       .groupby("gene")
#       .agg(
#           n_branches=("branch", "count"),
#           n_sig=("sig", "sum"),
#           sig_branches=("branch", lambda x: list(x[df.loc[x.index, "fdr"] < 0.05])),
#       )
# )

# summary["fraction_sig"] = summary["n_sig"] / summary["n_branches"]

# print(summary)
# summary.sort_values(by="fraction_sig", ascending=False).to_csv("absrel_summary_per_gene.csv")


# df_filtered = df[df["gene"].isin(sig_genes)]

# print("Гены, где все отмеченные ветви под отбором:")
# print(sig_genes)

# df.to_csv("absrel_full_table.csv", index=False)
# df_filtered.to_csv("absrel_genes_all_branches_selected.csv", index=False)

# print("Готово.")
