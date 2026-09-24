#!/usr/bin/env nextflow
/*
 * T2D Blood lncRNA Meta-Analysis Pipeline
 * =========================================
 * Standalone, T2D-only reanalysis extending Tkachenko et al. 2025 (IJMS,
 * DOI 10.3390/ijms262412046) with a genome-wide lncRNA layer and a
 * promoter TF-motif hypothesis-generation step (HOCOMOCO v11).
 *
 * Raw GEO data is downloaded separately (see README) into data/raw/ before
 * running this pipeline — GEO/EBI hosts sometimes require patient retry
 * logic (see download_cohorts.sh) that is easier to manage outside of
 * Nextflow's own retry semantics for large intermittent transfers.
 *
 * Usage:
 *   nextflow run main.nf
 */

nextflow.enable.dsl = 2

params.projectDir = "${workflow.projectDir}"
params.outdir     = "${params.projectDir}/data"

process FETCH_REFERENCE {
    publishDir "${params.outdir}/reference", mode: 'copy'
    input:
        path script
    output:
        path "gene_biotypes.csv"
        path "ensembl_to_symbol.csv"
    script:
    """
    cd ${params.projectDir}
    python3 ${script}
    cp data/reference/gene_biotypes.csv .
    cp data/reference/ensembl_to_symbol.csv .
    """
}

process PARSE_COHORTS {
    input:
        path script
        path biotype_ref
        path symbol_ref
    output:
        path "processed_marker.txt"
    script:
    """
    cd ${params.projectDir}
    python3 ${script}
    touch processed_marker.txt
    """
}

process META_ANALYSIS {
    publishDir "${params.outdir}/tables", mode: 'copy'
    input:
        path script
        path marker
    output:
        path "meta_analysis_deg_all_genes.csv"
        path "meta_analysis_deg_lncRNA.csv"
        path "meta_analysis_deg_protein_coding.csv"
        path "per_cohort_deg_all.csv"
    script:
    """
    cd ${params.projectDir}
    python3 ${script}
    cp data/tables/meta_analysis_deg_all_genes.csv .
    cp data/tables/meta_analysis_deg_lncRNA.csv .
    cp data/tables/meta_analysis_deg_protein_coding.csv .
    cp data/tables/per_cohort_deg_all.csv .
    """
}

process TOP_HITS_REPORT {
    input:
        path script
        path meta_all
        path meta_lnc
        path meta_pc
    output:
        path "report_done.txt"
    script:
    """
    cd ${params.projectDir}
    python3 ${script}
    touch report_done.txt
    """
}

process MOTIF_ANALYSIS {
    publishDir "${params.outdir}/tables", mode: 'copy'
    input:
        path script
        path report_marker
    output:
        path "promoter_motif_hits.csv"
    script:
    """
    cd ${params.projectDir}
    python3 ${script}
    cp data/tables/promoter_motif_hits.csv .
    """
}

workflow {
    ref_script    = file("${params.projectDir}/00_fetch_reference.py")
    parse_script  = file("${params.projectDir}/01_parse_cohorts.py")
    meta_script   = file("${params.projectDir}/02_meta_analysis_deg.py")
    report_script = file("${params.projectDir}/03_report_top_hits.py")
    motif_script  = file("${params.projectDir}/04_motif_analysis.py")

    (biotypes, symbols) = FETCH_REFERENCE(ref_script)
    marker = PARSE_COHORTS(parse_script, biotypes, symbols)
    (meta_all, meta_lnc, meta_pc, per_cohort) = META_ANALYSIS(meta_script, marker)
    report_marker = TOP_HITS_REPORT(report_script, meta_all, meta_lnc, meta_pc)
    MOTIF_ANALYSIS(motif_script, report_marker)
}
