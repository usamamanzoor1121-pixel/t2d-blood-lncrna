#!/usr/bin/env python3
"""
03_report_top_hits.py
=======================
Reports top nominal meta-analysis hits (by uncorrected p-value / concordance-
weighted effect score) regardless of whether they survive FDR correction —
useful for hypothesis generation when (as observed) genome-wide-significant
hits are rare/absent in a heterogeneous cross-study meta-analysis.

Also checks the specific genes from Wang et al. 2026 (PLOS One) BIRTHS-adjacent
T2D lncRNA panel and the BRS/BIRTHS candidate gene sets used elsewhere in this
project, for direct visibility into what this meta-analysis does/doesn't show
for those literature-nominated genes.
"""
import logging
from pathlib import Path
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

OUT = Path("data/tables")

# Wang et al. 2026 (PLOS One) 3-lncRNA panel gene/transcript IDs.
# NOTE: these are Ensembl transcript IDs (ENST) or StringTie novel-transcript IDs
# (MSTRG), not gene symbols. Our meta-analysis is gene-symbol-level (collapsed
# from transcript/Ensembl-gene level during parsing), so a direct ID match is not
# possible except in GSE280402, which is the only cohort we have at transcript
# resolution (ENST-indexed) BEFORE gene-level collapsing.
WANG_PANEL = ["ENST00000473095", "MSTRG.90147.1", "ENST00000531992"]

CANDIDATE_GENES = {
    "BIRTHS_stress": ["S100A8", "S100A9", "TXNIP", "IL1B", "NLRP3", "CCL2"],
    "BIRTHS_protect": ["IRS1", "PPARGC1A", "ADIPOR1", "IL10", "FOXO1", "TFAM"],
    "Tkachenko_consistent_down": ["FBLN2", "TPCN1", "PC", "SHANK1", "PLD4"],
}


def report_top_hits(meta_path: str, n=30, label=""):
    if not Path(meta_path).exists():
        log.warning(f"{meta_path} not found — skipping")
        return
    df = pd.read_csv(meta_path, index_col=0)
    log.info(f"\n{'='*70}\nTop {n} genes by meta p-value {label}\n{'='*70}")
    top = df.sort_values("pvalue_meta").head(n)
    cols = [c for c in ["n_studies", "mean_log2FC", "concordant_frac",
                         "pvalue_meta", "padj_meta", "cohorts"] if c in top.columns]
    log.info(f"\n{top[cols].to_string()}")
    return df


def check_candidates(df: pd.DataFrame, label=""):
    log.info(f"\n{'='*70}\nLiterature-candidate gene status in this meta-analysis {label}\n{'='*70}")
    for group, genes in CANDIDATE_GENES.items():
        log.info(f"\n[{group}]")
        present = df.loc[df.index.intersection(genes)]
        if present.empty:
            log.info("  none of these genes reached >=3 cohorts for meta-analysis")
        else:
            cols = [c for c in ["n_studies", "mean_log2FC", "concordant_frac",
                                 "pvalue_meta", "padj_meta"] if c in present.columns]
            log.info(f"\n{present[cols].to_string()}")


def check_wang_panel():
    log.info(f"\n{'='*70}\nWang et al. 2026 3-lncRNA panel — replication check\n{'='*70}")
    per_cohort = Path("data/tables/per_cohort_deg_all.csv")
    if not per_cohort.exists():
        log.warning("per_cohort_deg_all.csv not found")
        return
    deg = pd.read_csv(per_cohort)
    gse280402 = deg[deg["cohort"] == "GSE280402"]
    log.info(f"GSE280402 (only cohort with transcript-level ENST resolution "
             f"BEFORE gene-symbol collapsing): {len(gse280402)} genes tested")
    log.info("NOTE: our processed GSE280402_expr.csv is already collapsed to gene "
             "SYMBOL (summed across transcripts) during parsing — the specific "
             "Wang et al. transcript IDs (ENST00000473095, MSTRG.90147.1, "
             "ENST00000531992) are novel/specific transcript identifiers that "
             "cannot be recovered from symbol-level data. A true replication "
             "check would require re-parsing GSE280402 at transcript level and "
             "either (a) matching ENST00000473095/ENST00000531992 directly by ID "
             "[GSE280402 does have raw ENST IDs before collapsing], or (b) noting "
             "MSTRG.90147.1 is a StringTie-assembled novel transcript specific to "
             "Wang et al.'s own assembly and has no equivalent ID in GENCODE-"
             "annotated public data at all — it cannot be checked in ANY external "
             "cohort without re-running their exact assembly pipeline.")

    # Check what we CAN check: the two real Ensembl transcript IDs directly in
    # GSE280402's raw (pre-collapse) transcript-level file.
    raw = Path("data/raw/GSE280402/counts.tsv.gz")
    if raw.exists():
        raw_df = pd.read_csv(raw, sep="\t", compression="gzip", low_memory=False)
        for tid in ["ENST00000473095", "ENST00000531992"]:
            match = raw_df[raw_df["ENSEMBL"].str.startswith(tid)]
            if not match.empty:
                log.info(f"\n  {tid} found in GSE280402 raw transcript data:")
                log.info(f"  {match.to_string(index=False)}")
            else:
                log.info(f"\n  {tid} NOT present in GSE280402 (transcript not "
                         f"detected/quantified in this cohort)")


if __name__ == "__main__":
    df_all = report_top_hits(OUT / "meta_analysis_deg_all_genes.csv", label="(all biotypes)")
    if df_all is not None:
        check_candidates(df_all, label="(all-biotype meta-analysis)")

    lnc_path = OUT / "meta_analysis_deg_lncRNA.csv"
    if lnc_path.exists():
        report_top_hits(lnc_path, label="(lncRNA only)")

    check_wang_panel()
