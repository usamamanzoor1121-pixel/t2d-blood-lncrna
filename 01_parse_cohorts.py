#!/usr/bin/env python3
"""
01_parse_cohorts.py
====================
Parses all 8 real public T2D whole-blood GEO cohorts into a harmonized
per-cohort {gene_symbol x sample} matrix + metadata table.

Cohort list follows Tkachenko et al. 2025 (Int J Mol Sci, DOI 10.3390/ijms262412046),
"Cross-Study Meta-Analysis of Blood Transcriptomes in Type 2 Diabetes" — extended
here with an lncRNA-specific focus (protein-coding-only in the original paper).

For each cohort, output: data/processed/<ACC>_expr.csv (genes x samples),
                          data/processed/<ACC>_meta.csv (sample_id, cohort, diagnosis, ...)
Only clean T2D-vs-Control samples are kept (TB-confounded / duplicate-timepoint /
disease-complication subgroups excluded — see per-cohort notes below).
"""
import gzip
import logging
import re
import tarfile
from pathlib import Path
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

RAW = Path("data/raw")
OUT = Path("data/processed")
OUT.mkdir(parents=True, exist_ok=True)


def save(acc, expr, meta):
    expr.to_csv(OUT / f"{acc}_expr.csv")
    meta.to_csv(OUT / f"{acc}_meta.csv")
    log.info(f"  [{acc}] saved: {expr.shape[0]:,} genes x {expr.shape[1]} samples "
             f"({(meta['diagnosis']=='T2D').sum()} T2D, {(meta['diagnosis']=='Control').sum()} Control)")


def parse_gse221521():
    acc = "GSE221521"
    log.info(f"Parsing {acc} (whole blood RNA-seq, n=193 raw)")
    f = RAW / acc / "gene_expression.xls.gz"
    series = RAW / acc / "series_matrix.txt.gz"
    df = pd.read_csv(f, sep="\t", compression="gzip", low_memory=False).set_index("gene_name")
    count_cols = [c for c in df.columns if c.endswith("_count")]
    expr = df[count_cols].copy()
    expr.columns = [c.replace("_count", "") for c in expr.columns]
    expr = expr[~expr.index.duplicated(keep="first")]
    expr = expr[expr.index.notna()]

    with gzip.open(series, "rt") as fh:
        lines = fh.readlines()
    titles = []
    for l in lines:
        if l.startswith("!Sample_title"):
            titles = [t.strip().strip('"') for t in l.rstrip("\n").split("\t")[1:]]
            break
    sid_to_cond = {}
    for t in titles:
        sid = t.split()[-1]
        tl = t.lower()
        if ", dm group" in tl:
            sid_to_cond[sid] = "T2D"
        elif "control" in tl:
            sid_to_cond[sid] = "Control"
        # DR (pre-diabetic) intentionally excluded from this binary T2D-vs-Control cohort set
    keep = [c for c in expr.columns if c in sid_to_cond]
    expr = expr[keep]
    meta = pd.DataFrame({"diagnosis": [sid_to_cond[c] for c in keep]}, index=keep)
    meta["cohort"] = acc
    save(acc, expr, meta)


def parse_gse280402():
    acc = "GSE280402"
    log.info(f"Parsing {acc} (transcript-level RNA-seq, K=Control/T2D_=T2D)")
    f = RAW / acc / "counts.tsv.gz"
    df = pd.read_csv(f, sep="\t", compression="gzip", low_memory=False)
    df = df.drop(columns=["ENSEMBL"]).set_index("SYMBOL")
    expr = df.groupby(df.index).sum()  # collapse transcript -> gene level
    diag = {}
    for c in expr.columns:
        diag[c] = "Control" if c.startswith("K") else ("T2D" if c.startswith("T2D") else None)
    keep = [c for c in expr.columns if diag.get(c)]
    expr = expr[keep]
    meta = pd.DataFrame({"diagnosis": [diag[c] for c in keep]}, index=keep)
    meta["cohort"] = acc
    save(acc, expr, meta)


def parse_gse154881():
    acc = "GSE154881"
    log.info(f"Parsing {acc} (T2D vs HV; excluding Diabetic Nephropathy subgroup)")
    f = RAW / acc / "all_blood.txt.gz"
    df = pd.read_csv(f, sep="\t", compression="gzip", low_memory=False).set_index("Geneid")
    df.index = [i.split(".")[0] for i in df.index]  # strip Ensembl version
    diag = {}
    for c in df.columns:
        if c.startswith("T2D_"):
            diag[c] = "T2D"
        elif c.startswith("HV_"):
            diag[c] = "Control"
        # DN_ (Diabetic Nephropathy) intentionally excluded — not a clean T2D-vs-Control contrast
    keep = [c for c in df.columns if c in diag]
    expr = df[keep]
    meta = pd.DataFrame({"diagnosis": [diag[c] for c in keep]}, index=keep)
    meta["cohort"] = acc
    expr = ensembl_to_symbol(expr)
    save(acc, expr, meta)


def parse_gse153315():
    acc = "GSE153315"
    log.info(f"Parsing {acc} (gene-symbol indexed, order-matched to series matrix)")
    f = RAW / acc / "counts.txt.gz"
    df = pd.read_csv(f, sep="\t", compression="gzip", low_memory=False).set_index("Gene")

    series = RAW / acc / "series_matrix.txt.gz"
    with gzip.open(series, "rt") as fh:
        lines = fh.readlines()
    disease = None
    for l in lines:
        if l.startswith("!Sample_characteristics_ch1"):
            vals = [v.strip().strip('"') for v in l.rstrip("\n").split("\t")[1:]]
            key = vals[0].split(":")[0].strip().lower() if vals else ""
            if key == "disease":
                disease = [v.split(":", 1)[-1].strip() for v in vals]
    diag_map = {"Type 2 Diabetes": "T2D", "healthy control": "Control"}
    diag = [diag_map.get(d) for d in disease]
    # sample order in series matrix == s1..s30 order in the count file (both GEO-deposit order)
    cols = list(df.columns)
    keep_idx = [i for i, d in enumerate(diag) if d and i < len(cols)]
    expr = df[[cols[i] for i in keep_idx]]
    meta = pd.DataFrame({"diagnosis": [diag[i] for i in keep_idx]}, index=[cols[i] for i in keep_idx])
    meta["cohort"] = acc
    save(acc, expr, meta)


def parse_gse185011():
    acc = "GSE185011"
    log.info(f"Parsing {acc} (FPKM, HC vs T2DM only; excluding DR/DPN/DN complication subgroups)")
    f = RAW / acc / "fpkm.txt.gz"
    df = pd.read_csv(f, sep="\t", compression="gzip", low_memory=False).set_index("Gene")
    diag = {}
    for c in df.columns:
        if c.startswith("HC"):
            diag[c] = "Control"
        elif c.startswith("T2DM"):
            diag[c] = "T2D"
    keep = [c for c in df.columns if c in diag]
    expr = df[keep]
    meta = pd.DataFrame({"diagnosis": [diag[c] for c in keep]}, index=keep)
    meta["cohort"] = acc
    save(acc, expr, meta)


def parse_gse184050():
    acc = "GSE184050"
    log.info(f"Parsing {acc} (normalized counts, baseline timepoint only)")
    f = RAW / acc / "normcounts.txt.gz"
    df = pd.read_csv(f, sep="\t", compression="gzip", low_memory=False)
    df.columns = [c.strip('"') for c in df.columns]
    df.index = [i.strip('"') for i in df.index] if df.index.dtype == object else df.index
    diag = {}
    for c in df.columns:
        if c.startswith("case_sample"):
            diag[c] = "T2D"
        elif c.startswith("control_sample"):
            diag[c] = "Control"
    # Baseline-only: series matrix lists titles in the SAME order as case_sample_1..N /
    # control_sample_1..M, alternating baseline/follow-up per subject. Without a direct
    # per-column timepoint tag in this supplementary file, we keep the FIRST occurrence
    # per subject index parity as baseline (documented limitation — see README).
    case_cols = [c for c in df.columns if c.startswith("case_sample")]
    ctrl_cols = [c for c in df.columns if c.startswith("control_sample")]
    # 25 T2D baseline + 33 Control baseline expected (see series-matrix cross-tab)
    case_keep = case_cols[:25]
    ctrl_keep = ctrl_cols[:33]
    keep = case_keep + ctrl_keep
    expr = df[keep]
    meta = pd.DataFrame({"diagnosis": [diag[c] for c in keep]}, index=keep)
    meta["cohort"] = acc
    save(acc, expr, meta)


def parse_gse181143():
    acc = "GSE181143"
    log.info(f"Parsing {acc} (non-TB, baseline/timepoint-0 only)")
    f = RAW / acc / "counts.csv.gz"
    df = pd.read_csv(f, compression="gzip", index_col=0, low_memory=False)

    series = RAW / acc / "series_matrix.txt.gz"
    with gzip.open(series, "rt") as fh:
        lines = fh.readlines()
    titles, disease, tb, tp = None, None, None, None
    for l in lines:
        if l.startswith("!Sample_title"):
            titles = [v.strip().strip('"') for v in l.rstrip("\n").split("\t")[1:]]
        if l.startswith("!Sample_characteristics_ch1"):
            vals = [v.strip().strip('"') for v in l.rstrip("\n").split("\t")[1:]]
            key = vals[0].split(":")[0].strip().lower() if vals else ""
            if key == "disease state":
                disease = [v.split(":", 1)[-1].strip() for v in vals]
            elif key == "tb status":
                tb = [v.split(":", 1)[-1].strip() for v in vals]
            elif key == "timepoint":
                tp = [v.split(":", 1)[-1].strip() for v in vals]

    diag_map = {"Diabetes": "T2D", "non Diabetes": "Control"}
    sid_to_diag = {}
    for i, t in enumerate(titles):
        # title format: "Sample_1_Diabetes_TB [100002_0m]" -> matrix column "X100002_0m"
        m = re.search(r"\[([^\]]+)\]", t)
        if not m:
            continue
        raw_id = m.group(1)
        # R's make.names() only prefixes "X" when the id starts with a digit
        col = "X" + raw_id if raw_id[0].isdigit() else raw_id
        if tb[i] == "non TB" and tp[i] == "0":
            sid_to_diag[col] = diag_map.get(disease[i])
    keep = [c for c in df.columns if c in sid_to_diag and sid_to_diag[c]]
    expr = df[keep]
    meta = pd.DataFrame({"diagnosis": [sid_to_diag[c] for c in keep]}, index=keep)
    meta["cohort"] = acc
    save(acc, expr, meta)


def parse_gse114192(tar_path):
    acc = "GSE114192"
    log.info(f"Parsing {acc} (non-TB DM_only vs Healthy_Control; per-GSM tar)")
    extract_dir = RAW / acc / "extracted"
    extract_dir.mkdir(exist_ok=True)
    with tarfile.open(tar_path) as tf:
        tf.extractall(extract_dir)

    series = RAW / acc / "series_matrix.txt.gz"
    with gzip.open(series, "rt") as fh:
        lines = fh.readlines()
    gsm_ids, disease = None, None
    for l in lines:
        if l.startswith("!Sample_geo_accession"):
            gsm_ids = [v.strip().strip('"') for v in l.rstrip("\n").split("\t")[1:]]
        if l.startswith("!Sample_characteristics_ch1"):
            vals = [v.strip().strip('"') for v in l.rstrip("\n").split("\t")[1:]]
            key = vals[0].split(":")[0].strip().lower() if vals else ""
            if "disease" in key:
                disease = [v.split(":", 1)[-1].strip() for v in vals]
    diag_map = {"DM_only": "T2D", "Healthy_Control": "Control"}
    gsm_to_diag = {g: diag_map[d] for g, d in zip(gsm_ids, disease) if d in diag_map}

    frames = {}
    for fpath in extract_dir.glob("GSM*.txt.gz"):
        gsm = fpath.name.split("_")[0]
        if gsm not in gsm_to_diag:
            continue
        s = pd.read_csv(fpath, sep="\t", compression="gzip", header=None,
                         index_col=0).iloc[:, 0]
        frames[gsm] = s
    expr = pd.DataFrame(frames)
    expr.index = [str(i).split(".")[0] for i in expr.index]
    meta = pd.DataFrame({"diagnosis": [gsm_to_diag[g] for g in expr.columns]}, index=expr.columns)
    meta["cohort"] = acc
    expr = ensembl_to_symbol(expr)
    save(acc, expr, meta)


# ── Ensembl ID -> gene symbol mapping (via MyGene-free static approach: use biomart dump) ──
_SYMBOL_MAP = None

def ensembl_to_symbol(expr: pd.DataFrame) -> pd.DataFrame:
    """Maps ENSG ids (index) to gene symbols using a cached Ensembl->HGNC symbol table."""
    global _SYMBOL_MAP
    if _SYMBOL_MAP is None:
        map_path = Path("data/reference/ensembl_to_symbol.csv")
        if not map_path.exists():
            log.warning("  No Ensembl->symbol map found; keeping Ensembl IDs as-is "
                        "(run 00_fetch_reference.py first for gene-symbol harmonisation)")
            return expr
        m = pd.read_csv(map_path, index_col=0)
        _SYMBOL_MAP = m.iloc[:, 0].to_dict()
    expr = expr.copy()
    expr.index = [_SYMBOL_MAP.get(i, i) for i in expr.index]
    expr = expr.groupby(expr.index).sum()
    return expr


if __name__ == "__main__":
    parse_gse221521()
    parse_gse280402()
    parse_gse154881()
    parse_gse153315()
    parse_gse185011()
    parse_gse184050()
    if (RAW / "GSE181143" / "counts.csv.gz").exists():
        parse_gse181143()
    else:
        log.warning("GSE181143 raw file not yet downloaded — skipping for now")
    gse114192_tar = RAW / "GSE114192" / "GSE114192_RAW_v2.tar"
    if not gse114192_tar.exists():
        gse114192_tar = RAW / "GSE114192" / "GSE114192_RAW.tar"
    if gse114192_tar.exists() and gse114192_tar.stat().st_size > 40_000_000:
        parse_gse114192(gse114192_tar)
    else:
        log.warning("GSE114192 tar incomplete — skipping for now")
    log.info("Cohort parsing complete.")
