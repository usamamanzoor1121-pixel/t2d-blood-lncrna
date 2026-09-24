#!/usr/bin/env python3
"""
02_meta_analysis_deg.py
========================
Per-cohort differential expression (T2D vs Control) followed by a cross-study
meta-analysis, following the method of Tkachenko et al. 2025 (Int J Mol Sci,
DOI 10.3390/ijms262412046) — inverse-normal (Stouffer's) p-value combination,
weighted by sqrt(n), plus a concordance-aware effect-size score — but applied
separately to lncRNA-biotype and protein-coding genes (their study covered
protein-coding only).

Design choice vs Tkachenko et al.: rather than ComBat-seq batch-correcting all
cohorts into one pooled matrix (which assumes comparable library prep/platform
across very heterogeneous studies — RNA-seq vs different pipelines, whole blood
vs PBMC, different countries), we run DE PER COHORT on its own appropriately
normalised data, then combine at the p-value/effect-size level. This avoids
forcing platform-incompatible data into a single normalisation and is the more
conservative, defensible choice for such a heterogeneous cohort set.
"""
import logging
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

PROC = Path("data/processed")
REF = Path("data/reference")
OUT = Path("data/tables")
OUT.mkdir(parents=True, exist_ok=True)

MIN_GROUP_N = 4  # minimum samples per arm to attempt DE in a cohort


def median_of_ratios_normalise(counts: pd.DataFrame) -> pd.DataFrame:
    log_counts = np.log(counts.replace(0, np.nan))
    log_geo_mean = log_counts.mean(axis=1)
    finite = log_geo_mean.replace([np.inf, -np.inf], np.nan).notna()
    ratios = log_counts.loc[finite].sub(log_geo_mean.loc[finite], axis=0)
    size_factors = np.exp(ratios.median(axis=0))
    size_factors = size_factors.replace(0, np.nan).fillna(1.0)
    return counts.div(size_factors, axis=1)


# Cohorts already on a normalised/relative scale (FPKM, DESeq2-normalized) —
# these get log2(x+1) only; raw-count cohorts get size-factor normalisation first.
ALREADY_NORMALISED = {"GSE185011", "GSE184050"}


def load_cohort(acc: str):
    expr_path = PROC / f"{acc}_expr.csv"
    meta_path = PROC / f"{acc}_meta.csv"
    if not expr_path.exists():
        return None
    expr = pd.read_csv(expr_path, index_col=0)
    meta = pd.read_csv(meta_path, index_col=0)
    expr = expr[expr.index.notna()]
    expr = expr.groupby(expr.index).sum()  # collapse any remaining duplicate symbols
    return expr, meta


def cohort_deg(acc: str, expr: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    case_idx = meta[meta["diagnosis"] == "T2D"].index.intersection(expr.columns)
    ctrl_idx = meta[meta["diagnosis"] == "Control"].index.intersection(expr.columns)
    if len(case_idx) < MIN_GROUP_N or len(ctrl_idx) < MIN_GROUP_N:
        log.warning(f"  [{acc}] too few samples per arm (T2D={len(case_idx)}, "
                    f"Control={len(ctrl_idx)}) — skipping DE")
        return None

    e = expr[list(case_idx) + list(ctrl_idx)]
    e = e.loc[(e.sum(axis=1) > 0)]

    if acc not in ALREADY_NORMALISED:
        e = median_of_ratios_normalise(e)
    e_log = np.log2(e + 1)

    case_vals = e_log[case_idx].values
    ctrl_vals = e_log[ctrl_idx].values
    results = []
    for i, gene in enumerate(e_log.index):
        cv, kv = case_vals[i], ctrl_vals[i]
        if cv.std() + kv.std() < 1e-9:
            continue
        t, p = stats.ttest_ind(cv, kv, equal_var=False)
        lfc = cv.mean() - kv.mean()
        results.append({"gene": gene, "log2FC": lfc, "pvalue": p, "t_stat": t})
    df = pd.DataFrame(results)
    if df.empty:
        return None
    df["cohort"] = acc
    df["n_case"] = len(case_idx)
    df["n_ctrl"] = len(ctrl_idx)
    log.info(f"  [{acc}] DE run: {len(df):,} genes tested, "
             f"n_T2D={len(case_idx)}, n_Control={len(ctrl_idx)}")
    return df


def inverse_normal_combine(group: pd.DataFrame) -> pd.Series:
    """Stouffer's method: combine per-cohort p-values into one meta z/p, weighted
    by sqrt(n_total) of each cohort, using the SIGN of each cohort's t-statistic."""
    n_studies = len(group)
    weights = np.sqrt(group["n_case"] + group["n_ctrl"])
    z = stats.norm.isf(group["pvalue"] / 2) * np.sign(group["t_stat"])
    z_meta = (weights * z).sum() / np.sqrt((weights ** 2).sum())
    p_meta = 2 * stats.norm.sf(abs(z_meta))

    # Concordance-aware mean effect size: average log2FC, but flag if signs disagree
    signs = np.sign(group["log2FC"])
    concordant_frac = max((signs > 0).mean(), (signs < 0).mean())
    weighted_lfc = (weights * group["log2FC"]).sum() / weights.sum()
    # penalise discordant studies (matches the spirit of Tkachenko et al.'s
    # -1..1 concordance score, simplified to a multiplicative discount)
    effect_score = weighted_lfc * (2 * concordant_frac - 1)

    return pd.Series({
        "n_studies": n_studies,
        "z_meta": z_meta,
        "pvalue_meta": p_meta,
        "mean_log2FC": weighted_lfc,
        "concordant_frac": concordant_frac,
        "effect_score": effect_score,
        "cohorts": ",".join(sorted(group["cohort"].unique())),
    })


def run_meta_analysis(min_studies: int = 3):
    cohorts = ["GSE221521", "GSE280402", "GSE154881", "GSE153315",
               "GSE185011", "GSE184050", "GSE181143", "GSE114192"]

    all_deg = []
    for acc in cohorts:
        loaded = load_cohort(acc)
        if loaded is None:
            log.warning(f"  [{acc}] not available yet — skipping")
            continue
        expr, meta = loaded
        deg = cohort_deg(acc, expr, meta)
        if deg is not None:
            all_deg.append(deg)

    if not all_deg:
        log.error("No cohort DE results available — aborting")
        return

    combined = pd.concat(all_deg, ignore_index=True)
    combined.to_csv(OUT / "per_cohort_deg_all.csv", index=False)
    log.info(f"Per-cohort DE combined: {len(combined):,} gene-cohort rows "
             f"across {combined['cohort'].nunique()} cohorts")

    gene_counts = combined.groupby("gene")["cohort"].nunique()
    eligible_genes = gene_counts[gene_counts >= min_studies].index
    log.info(f"Genes tested in >= {min_studies} cohorts: {len(eligible_genes):,}")

    meta_rows = []
    for gene, grp in combined[combined["gene"].isin(eligible_genes)].groupby("gene"):
        row = inverse_normal_combine(grp)
        row["gene"] = gene
        meta_rows.append(row)
    meta_df = pd.DataFrame(meta_rows).set_index("gene")
    _, padj, _, _ = multipletests(meta_df["pvalue_meta"], method="fdr_bh")
    meta_df["padj_meta"] = padj
    meta_df = meta_df.sort_values("padj_meta")
    meta_df.to_csv(OUT / "meta_analysis_deg_all_genes.csv")

    n_sig = (meta_df["padj_meta"] < 0.05).sum()
    log.info(f"Meta-analysis significant genes (padj<0.05, all biotypes): {n_sig:,} / {len(meta_df):,}")

    # Split by biotype (lncRNA vs protein-coding) using GENCODE reference
    biotype_path = REF / "gene_biotypes.csv"
    if biotype_path.exists():
        biotypes = pd.read_csv(biotype_path, index_col=0)["gene_biotype"]
        meta_df["biotype"] = meta_df.index.map(biotypes)

        lnc = meta_df[meta_df["biotype"] == "lncRNA"]
        pc = meta_df[meta_df["biotype"] == "protein_coding"]
        lnc.to_csv(OUT / "meta_analysis_deg_lncRNA.csv")
        pc.to_csv(OUT / "meta_analysis_deg_protein_coding.csv")

        log.info(f"lncRNA genes in meta-analysis: {len(lnc):,}, "
                 f"significant (padj<0.05): {(lnc['padj_meta']<0.05).sum():,}")
        log.info(f"Protein-coding genes in meta-analysis: {len(pc):,}, "
                 f"significant (padj<0.05): {(pc['padj_meta']<0.05).sum():,}")
    else:
        log.warning("No biotype reference found — run 00_fetch_reference.py first")

    return meta_df


if __name__ == "__main__":
    run_meta_analysis()
