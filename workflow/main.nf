#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Multi-Omics Pipeline for Pentatrichomonas hominis Detection
 * Production-ready skeleton. Point --input at a real SRA manifest
 * and provide the custom DB built by scripts/build_custom_db.py.
 */

params.input          = "data/samples/sra_manifest.tsv"
params.outdir         = "results"
params.genome         = "data/genomes/Phominis_Hs3NIH.fna"
params.custom_db      = "results/detection/custom_db"
params.cpus           = 4
params.memory         = "8.GB"

log.info """
========================================================
 P. hominis Multi-Omics Detection Pipeline
========================================================
 input      : ${params.input}
 outdir     : ${params.outdir}
 genome     : ${params.genome}
 custom_db  : ${params.custom_db}
========================================================
"""

process BUILD_DB {
    tag "build_custom_db"
    publishDir "${params.outdir}/detection/custom_db", mode: 'copy'

    input:
    path genome
    path markers

    output:
    path "phominis_custom_db.fasta", emit: combined
    path "phominis_markers_only.fasta", emit: markers
    path "database_provenance.json"
    path "README_DB.md"

    script:
    """
    python ${projectDir}/../scripts/build_custom_db.py \
        --genome ${genome} \
        --markers ${markers} \
        --outdir .
    """
}

process DETECT {
    tag "${sample_id}"
    publishDir "${params.outdir}/detection/${sample_id}", mode: 'copy'
    cpus params.cpus
    memory params.memory

    input:
    tuple val(sample_id), path(reads)
    path db_dir

    output:
    path "detection_result.json"
    path "detection_report.txt"

    script:
    """
    python ${projectDir}/../scripts/detect_phominis.py \
        --fastq ${reads} \
        --db ${db_dir} \
        --outdir .
    """
}

workflow {
    // Placeholder: in real runs, use fromSRA or local FASTQs
    // For now the pipeline documents the intended structure.
    println "Nextflow workflow skeleton ready."
    println "1. Build DB with scripts/build_custom_db.py"
    println "2. Download real FASTQs via sra-tools using data/samples/sra_manifest.tsv"
    println "3. Run detection per sample"
    println "4. Aggregate and feed into scripts/ml_analysis.py"
}
