#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path
import re
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import chi2
import sys

# SLURM passes batch index
batch_file = sys.argv[1]
batch_id = Path(batch_file).stem.split("_")[1]

species=sys.argv[1]

PATH_TO_DATA = "/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/"
ALIGN_DIR = "alignments"

TREE = f"{PATH_TO_DATA}/initial_data/trees_labeled/tree_labelled_fg_inno_${species}.nwk"
TEMPLATE_ALT = f"{PATH_TO_DATA}/src/pair_analyis/codeml/template_alt.ctl"
TEMPLATE_NULL = f"{PATH_TO_DATA}/src/pair_analyis/codeml/template_null.ctl"

OUTDIR = Path(f"results_partial_{species}")
OUTDIR.mkdir(exist_ok=True)

tmpl_alt = open(TEMPLATE_ALT).read()
tmpl_null = open(TEMPLATE_NULL).read()

ogs_paths = [Path(x.strip()) for x in open(batch_file).readlines()]

def run_codeml(aln, og, mode):
    mode_dir = OUTDIR / og / mode
    mode_dir.mkdir(parents=True, exist_ok=True)

    outfile = mode_dir / "mlc"
    template = tmpl_alt if mode == "alt" else tmpl_null

    ctl_text = (template
                .replace("SEQFILE", str(aln.resolve()))
                .replace("TREEFILE", str(Path(TREE).resolve()))
                .replace("OUTFILE", str(outfile.resolve())))

    print(ctl_text)
    ctl_file = mode_dir / "codeml.ctl"
    ctl_file.write_text(ctl_text)

    with open(mode_dir/"codeml.log", "w") as log:
        subprocess.run(["codeml", "codeml.ctl"],
                       cwd=mode_dir,
                       stdout=log,
                       stderr=subprocess.STDOUT)

def parse_lnL(mlc_file):
    text = Path(mlc_file).read_text()
    m = re.search(r"lnL\(ntime:\s*\d+\s+np:\s*\d+\):\s*([-\d.]+)", text)
    return float(m.group(1)) if m else None

def parse_BEB_sites(mlc_file):
    text = Path(mlc_file).read_text()
    if "Bayes Empirical Bayes" not in text:
        return []

    beb_block = text.split("Bayes Empirical Bayes")[1]
    sites = []
    for line in beb_block.splitlines():
        m = re.match(r"\s*(\d+)\s+([A-Z*])\s+([0-9.]+)", line)
        if m:
            site = int(m.group(1))
            aa = m.group(2)
            pp = float(m.group(3))
            if pp >= 0.95:
                sites.append((site, aa, pp))
    return sites


rows = []

for aln in ogs_paths:
    og = aln.stem.replace(".codon_aln", "")

    run_codeml(aln, og, "alt")
    run_codeml(aln, og, "null")

    alt_mlc = OUTDIR / og / "alt" / "mlc"
    null_mlc = OUTDIR / og / "null" / "mlc"

    lnL_alt = parse_lnL(alt_mlc)
    print(lnL_alt)
    lnL_null = parse_lnL(null_mlc)
    print(lnL_null)

    # if lnL_alt is None or lnL_null is None:
    #     continue
    if lnL_alt is None or lnL_null is None:
        print(f"WARNING: lnL not parsed for {og}")
        continue

    LRT = 2 * (lnL_alt - lnL_null)
    p = 0.5 * (1 - chi2.cdf(LRT, df=1))

    beb = parse_BEB_sites(alt_mlc)
    beb_str = ";".join([f"{s}:{aa}:{pp:.3f}" for s, aa, pp in beb])

    rows.append([og, lnL_alt, lnL_null, LRT, p, beb_str])

df = pd.DataFrame(rows, columns=[
    "OG", "lnL_alt", "lnL_null", "LRT", "p_raw", "BEB_sites"
])

df.to_csv(f"results_partial/result_{batch_id}.tsv", sep="\t", index=False)
print("Batch done.")
