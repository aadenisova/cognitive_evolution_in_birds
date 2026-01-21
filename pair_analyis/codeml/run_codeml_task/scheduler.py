#!/usr/bin/env python3
import subprocess, time, re
from pathlib import Path
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import chi2

# ls alignments/*.codon_aln.phy | sed 's/.codon_aln.phy//' > og_list.txt

MAX_PARALLEL = 100    # максимальное число одновременно работающих codeml
SLEEP = 5             # интервал проверки очереди

PATH = "/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/"
TREE = f"{PATH}/codeml/tree_labelled_fg_inno.nwk"
TEMPLATE_ALT = f"{PATH}/src/pair_analyis/codeml/template_alt.ctl"
TEMPLATE_NULL = f"{PATH}/src/pair_analyis/codeml/template_null.ctl"

ogs = sorted([f for f in Path("alignments").glob("*.codon_aln.phy")])

active = {}     # jobID → (OG, ALIGN)
finished = []    # список OG с готовыми результатами

def launch(og_file):
    og = og_file.stem.replace(".codon_aln", "")
    cmd = [
        "sbatch", "--parsable", "run_one_og.sh",
        og, str(og_file.resolve()), TREE,
        TEMPLATE_ALT, TEMPLATE_NULL
    ]
    jobid = subprocess.check_output(cmd).decode().strip()
    active[jobid] = (og, og_file)
    print(f"LAUNCHED {og} → job {jobid}")

def is_done(jobid):
    out = subprocess.getoutput(f"squeue -j {jobid}")
    return len(out.splitlines()) <= 1  # нет в очереди → завершён

def parse_lnL(file):
    txt = file.read_text()
    m = re.search(r"lnL\(.*\):\s*([-\d.]+)", txt)
    return float(m.group(1)) if m else None

def parse_beb(file):
    txt = file.read_text()
    if "Bayes Empirical Bayes" not in txt:
        return ""
    block = txt.split("Bayes Empirical Bayes")[1]
    out = []
    for line in block.splitlines():
        m = re.match(r"\s*(\d+)\s+([A-Z*])\s+([0-9.]+)", line)
        if m and float(m.group(3)) >= 0.95:
            out.append(f"{m.group(1)}:{m.group(2)}:{m.group(3)}")
    return ";".join(out)

# MAIN LOOP --------------------------------------------------------------

queue = list(ogs)

results = []

print(f"Total OGs: {len(queue)}")

while queue or active:
    # 1) submit new jobs
    while queue and len(active) < MAX_PARALLEL:
        ogfile = queue.pop(0)
        launch(ogfile)

    # 2) check running jobs
    done = [jid for jid in active if is_done(jid)]

    for jid in done:
        og, ogfile = active[jid]
        print(f"FINISHED job {jid} ({og})")

        # parse results
        alt = Path("results")/og/"alt"/"mlc"
        null = Path("results")/og/"null"/"mlc"

        if alt.exists() and null.exists():
            lnL_alt = parse_lnL(alt)
            lnL_null = parse_lnL(null)

            if lnL_alt and lnL_null:
                LRT = 2*(lnL_alt - lnL_null)
                p = 0.5*(1 - chi2.cdf(LRT,1))
                beb = parse_beb(alt)

                results.append([og, lnL_alt, lnL_null, LRT, p, beb])

        del active[jid]

    time.sleep(SLEEP)

# export table
df = pd.DataFrame(results, columns=["OG","lnL_alt","lnL_null","LRT","p_raw","BEB"])
df["p_fdr"] = multipletests(df["p_raw"], method="fdr_bh")[1]
df.to_csv("branch_site_results.tsv", sep="\t", index=False)

print("DONE!")
