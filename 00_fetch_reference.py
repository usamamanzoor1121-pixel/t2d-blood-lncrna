#!/usr/bin/env python3
"""
00_fetch_reference.py
======================
Downloads the GENCODE human gene annotation (basic, comprehensive gene set)
and builds two reference tables used throughout the pipeline:
  data/reference/ensembl_to_symbol.csv   (ENSG id -> gene symbol, no version)
  data/reference/gene_biotypes.csv        (gene symbol -> GENCODE biotype)

The biotype table is what defines "lncRNA" downstream — using GENCODE's own
biotype calls (lncRNA, protein_coding, etc.) rather than an ad hoc heuristic.
"""
import gzip
import logging
import subprocess
from pathlib import Path
import re
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

REF = Path("data/reference")
REF.mkdir(parents=True, exist_ok=True)

GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
GTF_PATH = REF / "gencode.v44.annotation.gtf.gz"


def download():
    if GTF_PATH.exists() and GTF_PATH.stat().st_size > 30_000_000:
        log.info("GTF already downloaded")
        return
    log.info(f"Downloading GENCODE v44 GTF from {GTF_URL} (~45MB)...")
    subprocess.run(["curl", "-sS", "--retry", "8", "--retry-delay", "3",
                     "--retry-all-errors", "-o", str(GTF_PATH), GTF_URL], check=True)
    log.info(f"Downloaded: {GTF_PATH.stat().st_size/1e6:.1f} MB")


def parse_gtf():
    log.info("Parsing GTF gene records...")
    rows = []
    gene_id_re = re.compile(r'gene_id "([^"]+)"')
    gene_name_re = re.compile(r'gene_name "([^"]+)"')
    gene_type_re = re.compile(r'gene_type "([^"]+)"')

    with gzip.open(GTF_PATH, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            gid = gene_id_re.search(attrs)
            gname = gene_name_re.search(attrs)
            gtype = gene_type_re.search(attrs)
            if gid and gname and gtype:
                rows.append({
                    "gene_id": gid.group(1).split(".")[0],
                    "gene_symbol": gname.group(1),
                    "gene_biotype": gtype.group(1),
                })
    df = pd.DataFrame(rows).drop_duplicates(subset="gene_id")
    log.info(f"Parsed {len(df):,} gene records")
    log.info(f"Biotype counts:\n{df['gene_biotype'].value_counts().head(15)}")

    id_to_symbol = df.set_index("gene_id")["gene_symbol"]
    id_to_symbol.to_csv(REF / "ensembl_to_symbol.csv", header=["gene_symbol"])
    log.info(f"Saved: {REF / 'ensembl_to_symbol.csv'} ({len(id_to_symbol):,} entries)")

    symbol_to_biotype = df.drop_duplicates(subset="gene_symbol").set_index("gene_symbol")["gene_biotype"]
    symbol_to_biotype.to_csv(REF / "gene_biotypes.csv", header=["gene_biotype"])
    log.info(f"Saved: {REF / 'gene_biotypes.csv'} ({len(symbol_to_biotype):,} entries)")

    lncrna_biotypes = ["lncRNA"]
    n_lnc = (df["gene_biotype"].isin(lncrna_biotypes)).sum()
    log.info(f"GENCODE lncRNA-biotype genes: {n_lnc:,}")


if __name__ == "__main__":
    download()
    parse_gtf()
    log.info("Reference build complete.")
