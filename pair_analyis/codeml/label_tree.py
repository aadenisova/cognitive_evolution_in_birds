#!/usr/bin/env python3
"""
label_tree.py <tree_in.nwk> <comma_separated_taxon_list> <tree_out.nwk>

Adds #1 immediately after any matching taxon name (before ':') in the Newick string.
Matches exact token. If a taxon is absent, warns and continues.
"""
import sys, re
from pathlib import Path

if len(sys.argv) < 4:
    print("Usage: label_tree.py tree_in.nwk taxon1,taxon2,... tree_out.nwk")
    sys.exit(1)

tree_in = Path(sys.argv[1]).read_text().strip()
taxa = sys.argv[2].split(",")
out = sys.argv[3]

s = tree_in
missing = []
for t in taxa:
    # Replace occurrences like TaxonName: or TaxonName) or TaxonName, or TaxonName;
    # but only insert #1 before ':' if present, otherwise before ')' or ',' or ';'
    # We'll try replace "TaxonName:" -> "TaxonName#1:"; and "TaxonName)" -> "TaxonName#1)" etc.
    if t in s:
        # do multiple replacements safely
        s = re.sub(r'(?<![#\w])' + re.escape(t) + r'(?=[:\)\,;])', t + '{TEST}', s)
    else:
        missing.append(t)

if missing:
    print("Warning: the following taxa not found in tree:", missing)

Path(out).write_text(s + "\n")
print("Wrote labelled tree to", out)
