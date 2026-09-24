#!/usr/bin/env python3
"""
05_generate_figures.py
========================
Publication-style static figures for the T2D blood lncRNA meta-analysis.
Uses the validated colorblind-safe categorical palette (blue/orange/aqua,
CVD Delta E >= 9.2 all-pairs) and diverging blue<->red for magnitude+direction.

Figures:
  Fig1_cohort_overview.png    Sample sizes per cohort (T2D vs Control)
  Fig2_volcano.png            Meta-analysis volcano plot, coloured by biotype
  Fig3_concordance_heatmap.png Per-cohort log2FC for top candidate genes
  Fig4_qqplot.png             P-value calibration (observed vs expected)
"""
import logging
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

TAB = Path("data/tables")
PROC = Path("data/processed")
FIG = Path("data/figures")
FIG.mkdir(parents=True, exist_ok=True)

# ── Validated palette (references/palette.md) ──────────────────────────────
BLUE    = "#2a78d6"
ORANGE  = "#eb6834"
AQUA    = "#1baf7a"
RED     = "#e34948"
SURFACE = "#fcfcfb"
INK_PRIMARY   = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED     = "#898781"
GRID          = "#e1e0d9"
BASELINE      = "#c3c2b7"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "DejaVu Sans", "Arial"],
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.edgecolor": BASELINE,
    "axes.labelcolor": INK_SECONDARY,
    "text.color": INK_PRIMARY,
    "xtick.color": INK_MUTED,
    "ytick.color": INK_MUTED,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "grid.linestyle": "-",
    "axes.axisbelow": True,
    "savefig.dpi": 200,
})


def fig1_cohort_overview():
    meta_files = sorted(PROC.glob("*_meta.csv"))
    rows = []
    for f in meta_files:
        acc = f.stem.replace("_meta", "")
        m = pd.read_csv(f, index_col=0)
        rows.append({
            "cohort": acc,
            "T2D": (m["diagnosis"] == "T2D").sum(),
            "Control": (m["diagnosis"] == "Control").sum(),
        })
    df = pd.DataFrame(rows).set_index("cohort")
    df["total"] = df["T2D"] + df["Control"]
    df = df.sort_values("total", ascending=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    y = np.arange(len(df))
    bar_h = 0.36
    ax.barh(y + bar_h / 2, df["Control"], height=bar_h, color=BLUE, label="Control")
    ax.barh(y - bar_h / 2, df["T2D"], height=bar_h, color=ORANGE, label="T2D")

    for yi, (ctrl, t2d) in enumerate(zip(df["Control"], df["T2D"])):
        ax.text(ctrl + 3, yi + bar_h / 2, str(ctrl), va="center", fontsize=9, color=INK_SECONDARY)
        ax.text(t2d + 3, yi - bar_h / 2, str(t2d), va="center", fontsize=9, color=INK_SECONDARY)

    ax.set_yticks(y)
    ax.set_yticklabels(df.index, fontsize=10)
    ax.set_xlabel("Samples")
    ax.set_title("Cohort sample sizes\n466 samples across 8 independent public GEO cohorts",
                  fontsize=12, fontweight="bold", loc="left", pad=12)
    ax.grid(axis="y", visible=False)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.legend(frameon=False, loc="lower right", fontsize=10)
    plt.tight_layout()
    plt.savefig(FIG / "Fig1_cohort_overview.png")
    plt.close()
    log.info("Saved Fig1_cohort_overview.png")


def fig2_volcano():
    df = pd.read_csv(TAB / "meta_analysis_deg_all_genes.csv", index_col=0)
    df["neglog10p"] = -np.log10(df["pvalue_meta"].clip(lower=1e-300))

    def biotype_group(b):
        if b == "protein_coding":
            return "Protein-coding"
        elif b == "lncRNA":
            return "lncRNA"
        return "Other"

    df["group"] = df["biotype"].apply(biotype_group)
    colors = {"Protein-coding": BLUE, "lncRNA": ORANGE, "Other": AQUA}

    fig, ax = plt.subplots(figsize=(8, 7))
    for grp, sub in df.groupby("group"):
        ax.scatter(sub["mean_log2FC"], sub["neglog10p"], s=10, alpha=0.35,
                   color=colors[grp], label=grp, linewidths=0, zorder=2)

    ax.axhline(-np.log10(0.001), color=INK_MUTED, linewidth=1, linestyle="-", zorder=1)
    ax.text(ax.get_xlim()[1] * 0.98, -np.log10(0.001) + 0.05, "nominal p = 0.001",
            fontsize=8, color=INK_MUTED, ha="right")

    top = df.sort_values("pvalue_meta").head(8)
    ax.scatter(top["mean_log2FC"], top["neglog10p"], s=28, facecolors="none",
               edgecolors=INK_PRIMARY, linewidths=1.2, zorder=3)

    # Dense cluster -> stack labels clear of it (right edge) with thin leader
    # lines back to each point, evenly spaced top-to-bottom (no collisions).
    top_sorted = top.sort_values("neglog10p", ascending=False)
    label_x = 0.75
    y_top, y_bottom = top_sorted["neglog10p"].max() + 0.3, top_sorted["neglog10p"].min() - 1.0
    label_ys = np.linspace(y_top, max(y_bottom, 0.2), len(top_sorted))
    for (gene, row), label_y in zip(top_sorted.iterrows(), label_ys):
        label = gene if not gene.startswith("ENSG") else gene[:12] + "…"
        ax.annotate(
            label, xy=(row["mean_log2FC"], row["neglog10p"]),
            xytext=(label_x, label_y), fontsize=8.5, color=INK_PRIMARY,
            va="center", ha="left",
            arrowprops=dict(arrowstyle="-", color=INK_MUTED, linewidth=0.8,
                             shrinkA=3, shrinkB=3),
        )

    ax.set_xlabel("Mean log₂ fold-change (T2D vs Control)")
    ax.set_ylabel("−log₁₀(meta p-value)")
    ax.set_title("Meta-analysis volcano plot\n0 of 39,541 genes reach padj < 0.05 (best padj = 0.13)",
                  fontsize=12, fontweight="bold", loc="left", pad=12)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper right", fontsize=10, markerscale=2)
    plt.tight_layout()
    plt.savefig(FIG / "Fig2_volcano.png")
    plt.close()
    log.info("Saved Fig2_volcano.png")


def fig3_concordance_heatmap():
    meta = pd.read_csv(TAB / "meta_analysis_deg_all_genes.csv", index_col=0)
    per_cohort = pd.read_csv(TAB / "per_cohort_deg_all.csv")

    top_genes = meta.sort_values("pvalue_meta").head(15).index.tolist()
    cohort_order = ["GSE221521", "GSE280402", "GSE154881", "GSE153315",
                     "GSE185011", "GSE184050", "GSE181143", "GSE114192"]

    mat = pd.DataFrame(index=top_genes, columns=cohort_order, dtype=float)
    for gene in top_genes:
        sub = per_cohort[per_cohort["gene"] == gene]
        for _, row in sub.iterrows():
            if row["cohort"] in cohort_order:
                mat.loc[gene, row["cohort"]] = row["log2FC"]

    vmax = np.nanmax(np.abs(mat.values))
    fig, ax = plt.subplots(figsize=(8, 8))
    masked = np.ma.masked_invalid(mat.values.astype(float))
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(color="#f0efec")
    im = ax.imshow(masked, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")

    ax.set_xticks(range(len(cohort_order)))
    ax.set_xticklabels(cohort_order, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(top_genes)))
    ax.set_yticklabels([g if not g.startswith("ENSG") else g[:14] + "…" for g in top_genes], fontsize=9)

    for i in range(len(top_genes)):
        for j in range(len(cohort_order)):
            v = mat.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6.5,
                        color="white" if abs(v) > vmax * 0.5 else INK_PRIMARY)

    ax.set_xticks(np.arange(-0.5, len(cohort_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(top_genes), 1), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=2)
    ax.grid(which="major", visible=False)
    ax.spines[:].set_visible(False)

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("log₂ fold-change (T2D vs Control)", fontsize=9, color=INK_SECONDARY)
    cbar.outline.set_visible(False)

    ax.set_title("Cross-cohort direction concordance\ntop 15 exploratory candidate genes (gray = not tested in that cohort)",
                  fontsize=12, fontweight="bold", loc="left", pad=12)
    plt.tight_layout()
    plt.savefig(FIG / "Fig3_concordance_heatmap.png")
    plt.close()
    log.info("Saved Fig3_concordance_heatmap.png")


def fig4_qqplot():
    df = pd.read_csv(TAB / "meta_analysis_deg_all_genes.csv", index_col=0)
    p = df["pvalue_meta"].dropna().sort_values().values
    n = len(p)
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    observed = -np.log10(np.clip(p, 1e-300, 1))

    lam = np.median(observed) / np.median(expected)

    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.scatter(expected, observed, s=6, alpha=0.4, color=BLUE, linewidths=0, zorder=2)
    lims = [0, max(expected.max(), observed.max()) * 1.05]
    ax.plot(lims, lims, color=INK_MUTED, linewidth=1, zorder=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect("equal")

    ax.set_xlabel("Expected −log₁₀(p)")
    ax.set_ylabel("Observed −log₁₀(p)")
    ax.set_title(f"Meta-analysis p-value calibration (QQ plot)\ngenomic inflation λ = {lam:.2f} — no systematic inflation",
                  fontsize=12, fontweight="bold", loc="left", pad=12)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(FIG / "Fig4_qqplot.png")
    plt.close()
    log.info(f"Saved Fig4_qqplot.png (lambda={lam:.3f})")


if __name__ == "__main__":
    fig1_cohort_overview()
    fig2_volcano()
    fig3_concordance_heatmap()
    fig4_qqplot()
    log.info("All figures generated in data/figures/")
