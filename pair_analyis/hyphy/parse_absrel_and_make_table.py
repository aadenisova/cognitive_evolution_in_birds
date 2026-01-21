import json
import glob
import pandas as pd
import os
from statsmodels.stats.multitest import multipletests

# files = glob.glob("absrel_results_correct_inno_10/*.json")
files = glob.glob("absrel_results_inno_12/*.json")

rows = []

for f in files:
    gene = os.path.basename(f).replace(".json", "")

    if os.path.getsize(f) == 0:
        continue

    with open(f) as fh:
        data = json.load(fh)

    tested = data["tested"]["0"]
    test_branches = [b for b, t in tested.items() if t == "test"]

    branch_data = data["branch attributes"]["0"]

    for branch in test_branches:
        pval = branch_data[branch]["Uncorrected P-value"]
        rows.append({
            "gene": gene,
            "branch": branch,
            "raw_p": pval
        })

df = pd.DataFrame(rows)

# ---- FDR внутри каждой ВЕТВИ (вида) по всем генам ----
df["fdr"] = (
    df.groupby("branch")["raw_p"]
      .transform(lambda p: multipletests(p, method="fdr_bh")[1])
)

# Теперь фильтр: у каких генов ВСЕ тестируемые ветви значимы
significant_genes = (
    df.groupby("gene")
      .apply(lambda x: (x["fdr"] < 0.05).all())
)

sig_genes = significant_genes[significant_genes].index.tolist()

# Summary
summary = (
    df.assign(sig=df["fdr"] < 0.05)
      .groupby("gene")
      .agg(
          n_branches=("branch", "count"),
          n_sig=("sig", "sum"),
          sig_branches=("branch", lambda x: list(x[df.loc[x.index, "fdr"] < 0.05])),
      )
)

summary["fraction_sig"] = summary["n_sig"] / summary["n_branches"]
summary.sort_values(by="fraction_sig", ascending=False).to_csv("busted_e_summary_per_gene_filter.csv")

df_filtered = df[df["gene"].isin(sig_genes)]

print("Гены, где все ветви под FDR < 0.05:")
print(sig_genes)

df.to_csv("busted_e_full_table_filter.csv", index=False)
df_filtered.to_csv("busted_e_genes_all_branches_selected_filter.csv", index=False)

print("Готово.")


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
