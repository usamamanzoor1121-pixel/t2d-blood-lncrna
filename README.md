# T2D Blood lncRNA Meta-Analysis

A standalone, **Type 2 Diabetes–only** reanalysis of public whole-blood transcriptomic
data, extending [Tkachenko et al. 2025](https://doi.org/10.3390/ijms262412046)
(*Int J Mol Sci*, DOI 10.3390/ijms262412046) — a protein-coding-only meta-analysis of
8 public T2D blood cohorts — with a genome-wide **long non-coding RNA (lncRNA)**
layer they did not cover, plus a promoter transcription-factor motif
hypothesis-generation step. Also includes a partial replication check of
[Wang et al. 2026](https://doi.org/10.1371/journal.pone.0345359) (*PLOS One*)'s
proposed 3-lncRNA blood biomarker panel.

This project is intentionally **separate** from the earlier `blood-multiomics-ad-t2d`
project (AD × T2D) — T2D only, no AD, built from scratch with real public data from
the start (not a correction of prior fabricated results).

## Results

See **[`RESULTS_SUMMARY.md`](RESULTS_SUMMARY.md)** for the full write-up. Headline:

- **466 real samples** (244 T2D, 222 Control) pooled across 8 independent public GEO
  cohorts spanning 4+ countries, after excluding TB-comorbidity confounds,
  disease-complication subgroups, and duplicate longitudinal timepoints.
- **0 genes reach genome-wide FDR significance** (padj<0.05) under a conservative,
  replication-focused (Stouffer's inverse-normal) meta-analysis — a real, honestly
  reported null result that matches Tkachenko et al.'s own observation of "extremely
  low overlap between DEGs identified in different studies," using a stricter test
  of that same observation.
- A short list of **exploratory candidate genes** (EPHA2, ZBED3, TPCN1, and lncRNAs
  GABPB1-AS1, SIGLEC17P, LINC02044, among others) show strong nominal significance
  and 75-100% directional concordance across all 8 cohorts — reported as
  hypothesis-generating candidates for a dedicated, adequately-powered follow-up
  study, not as validated biomarkers.
- Promoter motif scan (HOCOMOCO v11) flags **NR1H4 (FXR)** motifs independently at
  two top candidates (TPCN1, SIGLEC17P) — a coherent, biologically plausible pattern
  worth mechanistic follow-up, not proof of regulation.
- Partial replication check of Wang et al. 2026's 3-lncRNA panel: 2 of 3 transcript
  IDs traced to real genes (CHCHD2, GRINA) in an independent cohort; the 3rd
  (a StringTie novel-assembly ID) is structurally unreplicable in any external dataset.

## Repository structure

```
T2D_lncRNA_Standalone/
├── 00_fetch_reference.py       # Downloads GENCODE v44 GTF, builds gene->biotype
│                                #   and Ensembl-ID->symbol reference tables
├── 01_parse_cohorts.py         # Parses 8 raw GEO downloads into harmonised
│                                #   {gene_symbol x sample} matrices + metadata
├── 02_meta_analysis_deg.py     # Per-cohort DE (Welch's t-test) + Stouffer's
│                                #   inverse-normal cross-study meta-analysis
├── 03_report_top_hits.py       # Top-hit reporting, literature-candidate check,
│                                #   Wang et al. 2026 panel replication check
├── 04_motif_analysis.py        # HOCOMOCO v11 promoter TF-motif scan
├── 04b_motif_scan_remaining.py # (helper — retries motif scan for genes whose
│                                #   sequence fetch initially timed out)
├── main.nf / nextflow.config   # Pipeline orchestration (see Reproducibility note)
├── setup_env.sh                # WSL Python venv setup
├── RESULTS_SUMMARY.md          # Full results write-up
└── data/
    ├── raw/           # Downloaded GEO data (not tracked in git — see .gitignore)
    ├── processed/     # Per-cohort harmonised expression + metadata
    ├── reference/     # GENCODE biotype table, HOCOMOCO motifs, promoter sequences
    ├── tables/        # DEG results, meta-analysis output, motif hits
    └── figures/       # (reserved for future visualisation)
```

## Cohorts

| Accession | T2D | Control | Platform | Notes |
|---|---|---|---|---|
| GSE221521 | 74 | 50 | RNA-seq | |
| GSE280402 | 8 | 8 | RNA-seq | Transcript-level (ENST) |
| GSE154881 | 5 | 5 | RNA-seq | Excl. Diabetic Nephropathy subgroup |
| GSE153315 | 20 | 10 | RNA-seq | |
| GSE185011 | 5 | 5 | RNA-seq (FPKM) | Excl. DR/DPN/DN subgroups |
| GSE184050 | 25 | 33 | RNA-seq (normalized) | Baseline timepoint only |
| GSE181143 | 55 | 75 | RNA-seq | Non-TB, baseline only (of 560 total) |
| GSE114192 | 52 | 36 | RNA-seq | Non-TB DM_only vs Healthy_Control (4 countries) |

## Reproducing this analysis

### 1. Environment (WSL, D: drive)

```bash
# from WSL, project lives at /mnt/d/T2D_lncRNA_Standalone
bash setup_env.sh   # creates venv/ with pandas, numpy, scipy, scikit-learn, etc.
source venv/bin/activate
```

### 2. Download raw data

Raw GEO files are not tracked in git (see `.gitignore`). Download each cohort's
supplementary/matrix files from NCBI GEO (accessions above) into
`data/raw/<ACCESSION>/` — see the download URLs and exact filenames referenced at
the top of `01_parse_cohorts.py`'s per-cohort parse functions. GSE114192 additionally
needs its `GSE114192_RAW.tar` extracted (handled automatically by the parser). The
GENCODE v44 GTF and HOCOMOCO v11 motif database are fetched automatically by
`00_fetch_reference.py` and `04_motif_analysis.py` respectively.

### 3. Run the pipeline

Either directly:

```bash
python3 00_fetch_reference.py
python3 01_parse_cohorts.py
python3 02_meta_analysis_deg.py
python3 03_report_top_hits.py
python3 04_motif_analysis.py
```

Or via Nextflow — **validated end-to-end** (all 5 processes succeeded, exit code 0,
17m40s, 0.6 CPU-hours; see `data/pipeline_report.html` / `pipeline_timeline.html`
for the execution report generated by that run). Requires a JDK ≥17 (WSL ships
Java 11 by default on some images — a portable JRE works fine and needs no root):

```bash
# if java -version shows <17, grab a portable JRE (no sudo required):
curl -sS -o /tmp/jdk17.tar.gz -L \
  'https://api.adoptium.net/v3/binary/latest/17/ga/linux/x64/jre/hotspot/normal/eclipse?project=jdk'
mkdir -p ~/jdk17 && tar -xzf /tmp/jdk17.tar.gz -C ~/jdk17 --strip-components=1
export JAVA_HOME=~/jdk17 PATH=~/jdk17/bin:$PATH

curl -s https://get.nextflow.io | bash   # installs ./nextflow
./nextflow run main.nf
```

Note: `04_motif_analysis.py`'s promoter-sequence fetches occasionally hit
transient Ensembl REST timeouts (external API flakiness, not a pipeline bug) — if
a single Nextflow run comes up short a gene or two in `promoter_motif_hits.csv`,
rerun `04b_motif_scan_remaining.py` for the missing genes and it will backfill them.

## Method notes

- **Per-cohort DE, not pooled/batch-corrected regression.** Each cohort is
  normalised on its own terms (median-of-ratios size-factor normalisation for raw
  RNA-seq counts; log2 only for already-normalised FPKM/DESeq2 cohorts), then tested
  independently (Welch's t-test, T2D vs Control). This is deliberately more
  conservative than ComBat-seq-pooling platform-incompatible studies together.
- **Stouffer's inverse-normal method** combines per-cohort p-values (weighted by
  √n), with a concordance-discounted effect size — a gene must show consistent
  cross-study signal, not just aggregate statistical power from 1-2 large cohorts.
- **Biotype classification** uses the full GENCODE v44 annotation (62,700 gene
  records, 18,866 lncRNA-biotype genes) — not a heuristic name-based filter.
- **Promoter motif scanning** uses HOCOMOCO v11 (Kulakovskiy et al.), scoring against
  each motif's own maximum-achievable PWM score (≥85% threshold) across both
  strands of a -2000/+200 bp window around each gene's annotated TSS.

## Author

Usama Manzoor · usama.manzoor1121@gmail.com
