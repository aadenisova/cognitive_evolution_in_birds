import sys, subprocess
from pathlib import Path

OG, ALIGN, TREE, TEMPLATE_ALT, TEMPLATE_NULL = sys.argv[1:]

OUT = Path("results")/OG

def run(template, mode):
    md = OUT/mode
    md.mkdir(parents=True, exist_ok=True)
    ctl = md/"codeml.ctl"
    mlc = md/"mlc"

    txt = open(template).read() \
        .replace("SEQFILE", ALIGN) \
        .replace("TREEFILE", TREE) \
        .replace("OUTFILE", str(mlc))

    ctl.write_text(txt)
    subprocess.run(["codeml", "codeml.ctl"], cwd=md)

run(TEMPLATE_ALT, "alt")
run(TEMPLATE_NULL, "null")
