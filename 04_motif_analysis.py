#!/usr/bin/env python3
"""
04_motif_analysis.py
======================
Mechanistic/regulatory layer: scans the promoters of the top exploratory
candidate genes from the meta-analysis for known transcription-factor binding
motifs, using HOCOMOCO v11 (Kulakovskiy et al., the TF motif database
maintained by V. Makeev's group, a frequent Medvedeva co-author — chosen
deliberately as the field-standard resource for this kind of question, not
an arbitrary choice).

This does NOT claim to demonstrate a regulatory mechanism — it generates
motif-based hypotheses (which TFs *could* plausibly bind these promoters)
for follow-up (ChIP-seq / reporter assay / ATAC-seq), consistent with the
exploratory, non-overclaiming framing of the rest of this project.
"""
import json
import logging
import re
from pathlib import Path
import requests
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

REF = Path("data/reference")
OUT = Path("data/tables")
OUT.mkdir(parents=True, exist_ok=True)

PROMOTER_UPSTREAM = 2000
PROMOTER_DOWNSTREAM = 200


def load_gene_coords():
    with open(REF / "motif_gene_coords.json") as f:
        return json.load(f)


def fetch_promoter_sequence(chrom, start, end, strand, gene):
    """Fetch -2000/+200 around the TSS, respecting strand, via Ensembl REST."""
    if strand == 1:
        tss = start
        region_start = max(1, tss - PROMOTER_UPSTREAM)
        region_end = tss + PROMOTER_DOWNSTREAM
    else:
        tss = end
        region_start = max(1, tss - PROMOTER_DOWNSTREAM)
        region_end = tss + PROMOTER_UPSTREAM

    url = (f"https://rest.ensembl.org/sequence/region/homo_sapiens/"
           f"{chrom}:{region_start}-{region_end}:{strand}?content-type=text/x-fasta")
    for attempt in range(5):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                lines = r.text.strip().split("\n")
                seq = "".join(lines[1:]).upper()
                log.info(f"  [{gene}] fetched promoter: {len(seq)} bp")
                return seq
        except Exception as e:
            log.warning(f"  [{gene}] attempt {attempt+1} failed: {e}")
    log.error(f"  [{gene}] FAILED to fetch promoter sequence")
    return None


def parse_meme_motifs(meme_path):
    """Parse MEME-format PWMs into {name: probability_matrix (Lx4, ACGT)}."""
    motifs = {}
    with open(meme_path) as f:
        content = f.read()
    blocks = content.split("MOTIF ")[1:]
    for block in blocks:
        lines = block.strip().split("\n")
        name = lines[0].split()[0]
        matrix_lines = []
        in_matrix = False
        for l in lines[1:]:
            if l.startswith("letter-probability matrix"):
                in_matrix = True
                continue
            if in_matrix:
                if not l.strip() or not re.match(r"^[\d.\-eE\s]+$", l.strip()):
                    break
                matrix_lines.append([float(x) for x in l.split()])
        if matrix_lines:
            motifs[name] = np.array(matrix_lines)
    log.info(f"Parsed {len(motifs)} motifs from {meme_path.name}")
    return motifs


_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def scan_sequence(seq: str, pwm: np.ndarray, bg=0.25, pseudocount=0.01):
    """Slide a PWM (probability matrix, LxACGT) across seq, return max log-odds
    score and its position, scanning both strands."""
    L = pwm.shape[0]
    log_pwm = np.log2((pwm + pseudocount) / (1 + 4 * pseudocount) / bg)
    max_score = -np.inf
    max_pos, max_strand = None, None

    def score_at(s):
        idxs = [_BASE_IDX.get(b) for b in s]
        if any(i is None for i in idxs):
            return -np.inf
        return sum(log_pwm[i, idxs[i]] for i in range(L))

    comp = str.maketrans("ACGT", "TGCA")
    rc = seq.translate(comp)[::-1]

    for i in range(len(seq) - L + 1):
        sc = score_at(seq[i:i + L])
        if sc > max_score:
            max_score, max_pos, max_strand = sc, i, "+"
    for i in range(len(rc) - L + 1):
        sc = score_at(rc[i:i + L])
        if sc > max_score:
            max_score, max_pos, max_strand = sc, len(seq) - L - i, "-"

    max_possible = log_pwm.max(axis=1).sum()
    return max_score, max_score / max_possible if max_possible > 0 else 0, max_pos, max_strand


def run(score_threshold_frac=0.85, top_n_motifs_per_gene=10):
    coords = load_gene_coords()
    motifs = parse_meme_motifs(REF / "HOCOMOCOv11_meme.txt")

    results = []
    for gene, info in coords.items():
        seq = fetch_promoter_sequence(info["seq_region_name"], info["start"],
                                       info["end"], info["strand"], gene)
        if not seq:
            continue
        gene_hits = []
        for motif_name, pwm in motifs.items():
            raw, frac, pos, strand = scan_sequence(seq, pwm)
            if frac >= score_threshold_frac:
                gene_hits.append({
                    "gene": gene, "motif": motif_name, "tf": motif_name.split("_")[0],
                    "score_frac_of_max": round(frac, 3), "position": pos, "strand": strand,
                })
        gene_hits.sort(key=lambda x: -x["score_frac_of_max"])
        results.extend(gene_hits[:top_n_motifs_per_gene])
        log.info(f"  [{gene}] {len(gene_hits)} motifs above {score_threshold_frac:.0%} "
                 f"of max PWM score (top {min(len(gene_hits), top_n_motifs_per_gene)} kept)")

    import pandas as pd
    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values(["gene", "score_frac_of_max"], ascending=[True, False])
        df.to_csv(OUT / "promoter_motif_hits.csv", index=False)
        log.info(f"\nSaved {len(df)} motif hits across {df['gene'].nunique()} genes "
                 f"to promoter_motif_hits.csv")
        log.info(f"\n{df.to_string(index=False)}")
    else:
        log.warning("No motif hits found above threshold")
    return df


if __name__ == "__main__":
    run()
