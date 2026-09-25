#!/usr/bin/env nextflow
// DSA4262 Assignment 1 - Task 5: long read RNA-Seq pipeline
// minimap2 -> samtools -> QC -> bambu

nextflow.enable.dsl = 2

params.reads          = "$projectDir/fastq/*.fastq.gz"
params.genome         = "$projectDir/ref/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
params.gtf            = "$projectDir/ref/Homo_sapiens.GRCh38.91.gtf"
params.outdir         = "results"
params.no_annotation  = false     // Scenario 2: run bambu without annotations
params.ndr            = 0.5       // NDR used when no annotation is supplied
params.min_reads      = 100000    // QC threshold: minimum reads per sample
params.min_mapped_pct = 80        // QC threshold: minimum % of reads mapped

process MINIMAP2 {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(fastq)
    path genome

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:
    // direct RNA reads are strand-specific (-uf) and noisier, so a shorter seed (-k14)
    def preset = sample_id.contains('directRNA') ? '-ax splice -uf -k14' : '-ax splice'
    """
    minimap2 ${preset} -t ${task.cpus} ${genome} ${fastq} > ${sample_id}.sam
    """
}

process SAMTOOLS_SORT {
    tag "$sample_id"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.bam ${sam}
    samtools index ${sample_id}.bam
    """
}

process QC {
    tag "$sample_id"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}.qc.tsv"

    script:
    """
    samtools stats -@ ${task.cpus} ${bam} > ${sample_id}.stats
    total=\$(awk -F'\\t' '\$1=="SN" && \$2=="raw total sequences:" {print \$3}' ${sample_id}.stats)
    mapped=\$(awk -F'\\t' '\$1=="SN" && \$2=="reads mapped:" {print \$3}' ${sample_id}.stats)
    avglen=\$(awk -F'\\t' '\$1=="SN" && \$2=="average length:" {print \$3}' ${sample_id}.stats)
    maxlen=\$(awk -F'\\t' '\$1=="SN" && \$2=="maximum length:" {print \$3}' ${sample_id}.stats)
    pct=\$(awk -v m=\$mapped -v t=\$total 'BEGIN{ if (t>0) printf "%.2f", 100*m/t; else print "0" }')
    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" ${sample_id} \$total \$mapped \$pct \$avglen \$maxlen > ${sample_id}.qc.tsv
    """
}

process QC_SUMMARY {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path qc_files

    output:
    path "qc_summary.tsv"

    script:
    """
    printf "sample\\ttotal_reads\\tmapped_reads\\tmapped_pct\\tavg_len\\tmax_len\\tstatus\\n" > qc_summary.tsv
    cat ${qc_files} | awk -F'\\t' -v mr=${params.min_reads} -v mp=${params.min_mapped_pct} \\
        'BEGIN{OFS="\\t"} {status = (\$2 >= mr && \$4 >= mp) ? "PASS" : "FAIL"; print \$0, status}' \\
        | sort >> qc_summary.tsv
    cat qc_summary.tsv
    """
}

process BAMBU {
    publishDir "${params.outdir}/bambu", mode: 'copy'

    input:
    path bams
    path bais
    path genome
    path gtf

    output:
    path "extended_annotations.gtf", optional: true
    path "*.txt"
    path "se.rds"

    script:
    def anno = params.no_annotation ? 'NONE' : gtf
    """
    run_bambu.R ${genome} ${anno} ${params.ndr} ${task.cpus} ${bams}
    """
}

workflow {
    genome_ch = Channel.value(file(params.genome, checkIfExists: true))
    gtf_ch    = Channel.value(file(params.gtf,    checkIfExists: true))

    reads_ch = Channel
        .fromPath(params.reads, checkIfExists: true)
        .map { f -> tuple(f.simpleName, f) }

    sam_ch = MINIMAP2(reads_ch, genome_ch)
    bam_ch = SAMTOOLS_SORT(sam_ch)

    QC_SUMMARY(QC(bam_ch).collect())

    BAMBU(bam_ch.map { it[1] }.collect(),
          bam_ch.map { it[2] }.collect(),
          genome_ch,
          gtf_ch)
}
