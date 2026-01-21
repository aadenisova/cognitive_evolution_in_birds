#!/usr/bin/env python3
from pathlib import Path
import math
import os

ALIGN_DIR = Path("alignments")
MAX_ARRAY = 100  # ограничение SLURM

ogs = sorted([f for f in ALIGN_DIR.glob("*.macse.clean.phy") if f.stat().st_size > 0])
n = len(ogs)
batch_size = math.ceil(n / MAX_ARRAY)

Path("batches").mkdir(exist_ok=True)

for i in range(MAX_ARRAY):
    start = i * batch_size
    end = start + batch_size
    sub = ogs[start:end]
    if not sub:
        break

    with open(f"batches/batch_{i:03d}.txt", "w") as f:
        for p in sub:
            f.write(str(p) + "\n")

print("Created batches.")
