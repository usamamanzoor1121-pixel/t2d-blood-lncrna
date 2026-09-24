#!/usr/bin/env python3
"""Scans the 3 promoter sequences fetched separately via curl (GABPB1-AS1,
SIGLEC17P, LINC02044) and appends results to promoter_motif_hits.csv."""
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, ".")
import importlib.util
spec = importlib.util.spec_from_file_location("motif_analysis", "04_motif_analysis.py")
ma = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ma)

REF = Path("data/reference")
OUT = Path("data/tables")

motifs = ma.parse_meme_motifs(REF / "HOCOMOCOv11_meme.txt")

genes = ["GABPB1-AS1", "SIGLEC17P", "LINC02044"]
results = []
for gene in genes:
    fa_path = REF / f"promoter_{gene}.fa"
    lines = fa_path.read_text().strip().split("\n")
    seq = "".join(lines[1:]).upper()
    print(f"[{gene}] {len(seq)} bp")

    gene_hits = []
    for motif_name, pwm in motifs.items():
        raw, frac, pos, strand = ma.scan_sequence(seq, pwm)
        if frac >= 0.85:
            gene_hits.append({
                "gene": gene, "motif": motif_name, "tf": motif_name.split("_")[0],
                "score_frac_of_max": round(frac, 3), "position": pos, "strand": strand,
            })
    gene_hits.sort(key=lambda x: -x["score_frac_of_max"])
    results.extend(gene_hits[:10])
    print(f"  {len(gene_hits)} motifs above 85% threshold, top 10 kept")

new_df = pd.DataFrame(results)
existing = pd.read_csv(OUT / "promoter_motif_hits.csv")
combined = pd.concat([existing, new_df], ignore_index=True)
combined = combined.sort_values(["gene", "score_frac_of_max"], ascending=[True, False])
combined.to_csv(OUT / "promoter_motif_hits.csv", index=False)
print(f"\nCombined total: {len(combined)} hits across {combined['gene'].nunique()} genes")
print(new_df.to_string(index=False))
