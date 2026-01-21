from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import json
import sys

#########################################
# ФАЙЛЫ – измени под себя
#########################################

og = "OG0007618"
FEL_FILE = f"{og}.json"        # или .csv
CODON_ALIGNMENT = f"alignments/{og}.macse.clean.renamed.fa"
OUTPUT = "fel_processed.tsv"

#########################################
# 1. Загрузить FEL-результат
#########################################

def load_fel(path):
    if path.endswith(".json"):
        with open(path) as f:
            data = json.load(f)
        rows = data["MLE"]["content"]["0"] #["tables"]["By Site"]
        df = pd.DataFrame(rows[1:], columns=rows[0])
        df = df.rename(columns={"Site":"site"})
        df["site"] = df["site"].astype(int)
        return df
    else:
        # CSV / TSV
        return pd.read_csv(path)

fel = load_fel(FEL_FILE)

#########################################
# 2. Загрузка кодонного выравнивания
#########################################

records = list(SeqIO.parse(CODON_ALIGNMENT, "fasta"))
taxa = [r.id for r in records]

def translate_codon_alignment(seq):
    seq = str(seq)
    aa = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) < 3 or "-" in codon:
            aa.append("-")
        else:
            aa.append(str(Seq(codon).translate()))
    return aa

translated = {r.id: translate_codon_alignment(r.seq) for r in records}
alignment_length = len(next(iter(translated.values())))

#########################################
# 3. Построить таблицу с аминокислотами по сайтам
#########################################

rows = []

for _, row in fel.iterrows():
    site = int(row["site"]) - 1   # HyPhy использует 1-based
    
    record = {
        "site": site+1,
    }
    
    # копируем численные поля FEL:
    for col in row.index:
        if col != "site":
            record[col] = row[col]
    
    # добавляем аминокислоты таксонов
    for taxon in taxa:
        aa = translated[taxon][site]
        record[f"AA_{taxon}"] = aa
    
    # консенсус
    aas = [translated[t][site] for t in taxa if translated[t][site] != "-"]
    record["AA_consensus"] = max(set(aas), key=aas.count) if aas else "-"
    
    rows.append(record)

df_out = pd.DataFrame(rows)

#########################################
# 4. Сохранение
#########################################

df_out.to_csv(OUTPUT, sep="\t", index=False)

print(f"Готово! Результат сохранён в {OUTPUT}")
