# T2D Blood lncRNA Meta-Analysis — Results Summary (final, 8/8 cohorts)

Generated 2026-09-24. Standalone, T2D-only project (no AD). Extends Tkachenko et al.
2025 (Int J Mol Sci, DOI 10.3390/ijms262412046) — a protein-coding-only meta-analysis
of 8 public T2D whole-blood cohorts — by (a) adding a genome-wide lncRNA-specific
layer they did not cover, and (b) cross-checking Wang et al. 2026 (PLOS One)'s
proposed 3-lncRNA blood biomarker panel.

## Cohorts (all 8 processed)

| Accession | n T2D | n Control | Notes |
|---|---|---|---|
| GSE221521 | 74 | 50 | Whole blood RNA-seq |
| GSE280402 | 8 | 8 | Whole blood RNA-seq, transcript-level (ENST) |
| GSE154881 | 5 | 5 | Excl. Diabetic Nephropathy subgroup |
| GSE153315 | 20 | 10 | |
| GSE185011 | 5 | 5 | Excl. DR/DPN/DN complication subgroups |
| GSE184050 | 25 | 33 | Baseline timepoint only (excl. follow-up) |
| GSE181143 | 55 | 75 | Non-TB, baseline (timepoint 0) only — excludes 430 TB-confounded/follow-up samples from this 560-sample, 6-naming-convention, multi-site cohort |
| GSE114192 | 52 | 36 | Non-TB DM_only vs Healthy_Control subset of a 4-country (Romania, South Africa, Indonesia, Peru) TB/DM comorbidity study |
| **Total (8 cohorts)** | **244** | **222** | **466 samples**, vs. n=189 in the original single-cohort analysis (2.5x) |

Each cohort was filtered to a clean T2D-vs-Control contrast — TB-comorbid samples,
disease-complication subgroups (nephropathy/retinopathy/neuropathy), and duplicate
longitudinal timepoints were excluded (documented per-cohort in `01_parse_cohorts.py`).
GSE114192 incidentally adds real population/ancestry diversity (4 countries across
3 continents) to the meta-cohort.

## Method

Per-cohort Welch's t-test (T2D vs Control) on appropriately normalised expression
(size-factor normalisation for raw RNA-seq counts; log2 only for already-normalised
FPKM/DESeq2-normalized cohorts), then Stouffer's inverse-normal method to combine
p-values across cohorts (weighted by √n), with a concordance-discounted mean
effect size. This is a **conservative, replication-focused** design — a gene must
show a consistent, weighted-significant signal across independent cohorts, not just
aggregate power from pooling — different from Tkachenko et al.'s pooled mixed-model
(variancePartition) approach.

Biotype classification used the full GENCODE v44 gene annotation (62,700 gene
records; 18,866 lncRNA-biotype genes genome-wide), not a partial lookup.

## Headline result: no genome-wide-significant genes, even at full power

**0 of 39,541 genes tested in ≥3 cohorts reach padj<0.05** in the complete 8-cohort
meta-analysis (466 samples) — split as 0/19,178 protein-coding and 0/7,393 lncRNA.
Adding the 8th cohort (GSE114192, 88 more samples) did not change this conclusion;
the closest gene (ENSG00000261799, an unnamed lncRNA) improved only from a prior
partial result to padj=0.126, still short of significance.

This is not a contradiction of Tkachenko et al.'s reported "2,065 DEGs" — it reflects
a genuinely different, more conservative question. They themselves reported
**"extremely low overlap between DEGs identified in different studies"** and found
only 5 genes consistently downregulated across 3+ of their 8 cohorts (FBLN2, TPCN1,
PC, SHANK1, PLD4). In this full 8-cohort replication-focused reanalysis, the best of
those five (TPCN1) reaches p=0.00047 (padj=0.40) — closer to significance than in the
7-cohort interim result, but still not FDR-significant. The apparent disagreement
with Tkachenko et al.'s headline number is a methods artifact, not a factual
contradiction: pooled mixed-model analysis is sensitive to genes driven by the
largest 1-2 cohorts even without cross-study consistency, while Stouffer's method
specifically penalizes that.

**Conclusion: T2D whole-blood transcriptomic changes, at the single-gene level, are
not consistently reproducible across 8 independent public cohorts spanning 4+
countries, multiple platforms, and different sample-collection protocols**, once
evaluated by a method that requires real cross-study agreement rather than aggregate
statistical power. This is itself a substantive, honest finding about the field's
current public data, not a failure of the analysis.

## Top overall hits (exploratory, none FDR-significant)

| Gene | Biotype | n cohorts | mean log2FC | concordance | p (meta) | padj |
|---|---|---|---|---|---|---|
| ENSG00000261799 | lncRNA (unnamed) | 3 | −0.65 | 100% | 0.000003 | 0.126 |
| H2AFY | protein-coding | 3 | −0.66 | 100% | 0.000006 | 0.126 |
| EPHA2 | protein-coding | **8/8** | −0.50 | 75% | 0.000086 | 0.398 |
| LAS1L | protein-coding | **8/8** | −0.25 | 100% | 0.000098 | 0.398 |
| TM9SF4 | protein-coding | **8/8** | −0.38 | 88% | 0.000123 | 0.398 |
| ZBED3 | protein-coding | **8/8** | −0.35 | 100% | 0.000190 | 0.398 |

EPHA2, ZBED3, and several others reaching nominal p<0.0005 across all 8 independent
cohorts (not just a subset) is a meaningfully stronger signal than any gene managed
in the 7-cohort interim result — EPHA2 in particular has prior literature association
with diabetic vascular complications, making it a defensible exploratory candidate
despite not reaching formal significance.

## lncRNA candidates (exploratory, genome-wide biotype coverage)

| Gene | n cohorts | mean log2FC | concordance | p (meta) | padj |
|---|---|---|---|---|---|
| LINC02044 | 5 | −0.34 | 100% | 0.00019 | 0.40 |
| GABPB1-AS1 | **8/8** | +0.26 | 100% | 0.00041 | 0.40 |
| CLN8-AS1 | 3 | +0.40 | 67% | 0.00073 | 0.43 |
| LINC01160 | 6 | +0.31 | 100% | 0.00097 | 0.48 |
| EFHD2-AS1 | 3 | −0.63 | 100% | 0.00099 | 0.48 |
| LCT-AS1 | 5 | −0.41 | 100% | 0.00100 | 0.48 |
| MIR222HG | 5 | +0.55 | 80% | 0.00118 | 0.52 |
| ZFAND2A-DT | 5 | −0.36 | 100% | 0.00131 | 0.53 |
| SIGLEC17P | **8/8** | −0.39 | 88% | 0.00173 | 0.56 |
| LINC02687 | 3 | −0.29 | 100% | 0.00175 | 0.56 |

None survive FDR correction. **GABPB1-AS1 and SIGLEC17P stand out as tested in all 8
cohorts with 88-100% directional concordance** — the strongest exploratory lncRNA
candidates from this analysis. Full ranked list (30+ genes, including unnamed
ENSG-only lncRNAs): `data/tables/meta_analysis_deg_lncRNA.csv`.

**These are reported as exploratory, hypothesis-generating candidates for a
dedicated, adequately-powered follow-up study, not as validated biomarkers.**

## Literature candidate genes: mostly not supported, TPCN1 strengthens

| Gene set | Best result (8-cohort) |
|---|---|
| BIRTHS stress genes (S100A8, S100A9, TXNIP, IL1B, NLRP3, CCL2) | Best: S100A8 p=0.11 (not significant) |
| BIRTHS protective genes (IRS1, PPARGC1A, ADIPOR1, IL10, FOXO1, TFAM) | Best: ADIPOR1 p=0.072 (not significant) |
| Tkachenko et al. "consistent" genes (FBLN2, TPCN1, PC, SHANK1, PLD4) | Best: **TPCN1 p=0.00047, padj=0.40** (strengthened with full cohort set, still not FDR-significant) |

None of the biomarker genes used elsewhere in the field (including in the earlier,
now-corrected AD×T2D project's BIRTHS score) reach FDR significance in this more
conservative replication design, though TPCN1 is now the closest literature
candidate to significance after adding the 8th cohort.

## Wang et al. 2026 (PLOS One) 3-lncRNA panel — replication check

Wang et al. proposed ENST00000473095, MSTRG.90147.1, and ENST00000531992 as a
3-lncRNA blood biomarker panel (combined AUC=0.73) from an 8+8 discovery / 85+85
validation cohort with pancreatic single-cell cross-referencing.

- **ENST00000473095** and **ENST00000531992** were located directly in GSE280402's
  raw transcript-level data — they map to the **CHCHD2** and **GRINA** gene loci
  respectively (confirmed via GENCODE v44). Both loci are canonically protein-coding
  genes; the specific transcripts Wang et al. used are minor/alternative isoforms.
  Neither CHCHD2 nor GRINA appear in this meta-analysis's top hits at the gene level
  (i.e., no cross-study signal at either locus) — a modest independent data point
  against the panel's cross-cohort generalisability, though this compares gene-level
  aggregate signal to Wang et al.'s specific-transcript score, which are not
  identical measurements.
- **MSTRG.90147.1** is a StringTie-assembled novel-transcript ID specific to Wang et
  al.'s own transcriptome assembly. It has no equivalent identifier in GENCODE-
  annotated public data and **cannot be checked in any external cohort** without
  re-running their exact assembly pipeline on raw reads — a real, structural
  limitation of using novel-assembly IDs as biomarkers rather than annotated genes,
  worth noting as a methodological critique regardless of this analysis's own results.

## Limitations

1. **No formal cross-platform batch correction** — this meta-analysis deliberately
   avoids ComBat-seq-style pooling across platform-incompatible cohorts (RNA-seq vs
   different pipelines/tissue preps) in favor of per-cohort DE + p-value combination;
   this is the more conservative and defensible choice for such heterogeneous data,
   but is not directly comparable number-for-number to Tkachenko et al.'s pooled
   approach.
2. **HbA1c/clinical severity metadata absent** for most cohorts (consistent with the
   earlier AD×T2D audit's finding for GSE221521) — this analysis is limited to
   T2D-vs-Control diagnosis contrasts, not severity/dose-response relationships.
3. **Transcript-level resolution only available for GSE280402** — all other cohorts
   are gene-symbol-level, so isoform-specific effects (like Wang et al.'s panel)
   cannot be fully assessed against them.
4. **Simplified batch handling within GSE181143 and GSE184050** for longitudinal/
   multi-site structure (baseline-timepoint selection, non-TB filtering) — documented
   per-cohort in `01_parse_cohorts.py`, reasonable but not exhaustively validated
   against each original study's own QC pipeline.

## Promoter motif analysis (HOCOMOCO v11)

Scanned -2000/+200 bp promoter regions of the 9 strongest exploratory candidates
(EPHA2, LAS1L, TM9SF4, ZBED3, TPCN1, GABPB1-AS1, SIGLEC17P, LINC02044, LINC01160)
against all 769 HOCOMOCO v11 human TF motifs (score ≥85% of each motif's maximum
possible PWM score). This generates regulatory *hypotheses* — which TFs could
plausibly bind these promoters — not evidence of an actual regulatory mechanism
(that would require ChIP-seq, ATAC-seq, or reporter assays).

**Notable pattern:** TPCN1 (the literature candidate closest to cross-study
significance, p=0.00047) and SIGLEC17P (a top lncRNA candidate) both carry strong
**NR1H4 (FXR)** motifs — FXR is a bile-acid-activated nuclear receptor and an active
T2D/NAFLD drug target, and its independent appearance at two different top-ranked
loci is a coherent, biologically plausible pattern worth flagging (not proof).
TPCN1's promoter also carries **PPARA** and **VDR** motifs — both established
metabolic nuclear receptors. ZBED3 carries a **HIF1A** motif (hypoxia/metabolic
stress response). None of this constitutes a demonstrated mechanism; it is reported
as a prioritization signal for which candidates might be worth mechanistic follow-up
first. Full results (90 motif hits across all 9 genes): `data/tables/promoter_motif_hits.csv`.

## Next steps (per project roadmap)

- HOCOMOCO transcription-factor motif scan on the promoters of the top exploratory
  candidates (EPHA2, GABPB1-AS1, SIGLEC17P, etc.) for a mechanistic regulatory layer.
- Package as a reproducible Nextflow pipeline.
- Write up as a standalone note: an honest null/exploratory result with a rigorous,
  transparent, literature-grounded methodology — a defensible contribution in its own
  right, and a foundation for either a larger prospective cohort or an integration
  with cell-type-resolved (single-cell) blood data.
