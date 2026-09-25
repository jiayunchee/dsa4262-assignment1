# DSA4262 Assignment 1 — Long-read RNA-Seq pipeline

A Nextflow (DSL2) pipeline for Oxford Nanopore long-read RNA-Seq: alignment,
BAM conversion, quality control, and transcript discovery and quantification.
Built for DSA4262 Assignment 1 using four samples from the
[SG-NEx](https://github.com/GoekeLab/sg-nex-data) dataset.

## Pipeline steps

| Process | Tool | Description |
|---|---|---|
| `MINIMAP2` | minimap2 2.31 | Spliced alignment to the genome. Flags are chosen per sample: direct RNA uses `-ax splice -uf -k14`, cDNA uses `-ax splice`. |
| `SAMTOOLS_SORT` | samtools | Converts SAM to a coordinate-sorted BAM and indexes it. |
| `QC` | samtools stats | Per-sample read count, mapped reads, mapping rate, mean and maximum read length. |
| `QC_SUMMARY` | awk | Merges the per-sample metrics and flags each sample PASS or FAIL against the thresholds. |
| `BAMBU` | bambu 3.12.1 | Transcript discovery and quantification, with or without a reference annotation. |

**Inputs:** fastq (one per sample), genome fasta, annotation GTF
**Outputs:** sorted and indexed BAM files, extended annotation GTF, transcript
and gene count tables, QC summary

## Why the alignment flags differ

Direct RNA reads are sequenced from native RNA, so the transcript strand is
known: `-uf` restricts alignment to the forward transcript strand. Their higher
error rate is handled by a shorter seed length, `-k14`. cDNA reads can originate
from either strand, so neither flag applies. Both protocols need `-ax splice`,
because RNA reads span introns.

## Requirements

- Nextflow 24+
- Two conda environments, with their paths set in `nextflow.config`:
  - `dsa4262`: minimap2, samtools
  - `bambu`: R with the bambu Bioconductor package

```bash
conda create -n dsa4262 -c conda-forge -c bioconda minimap2 samtools nextflow
conda create -n bambu   -c conda-forge -c bioconda bioconductor-bambu bioconductor-rtracklayer r-biocmanager
```

## Running

Scenario 1 — Bambu with genome annotations:

```bash
nextflow run main.nf \
  --reads '/path/to/fastq/*.fastq.gz' \
  --genome /path/to/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  --gtf /path/to/Homo_sapiens.GRCh38.91.gtf \
  -with-report report_s1.html
```

Scenario 2 — Bambu without annotations, reusing cached alignments:

```bash
nextflow run main.nf \
  --reads '/path/to/fastq/*.fastq.gz' \
  --genome /path/to/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  --gtf /path/to/Homo_sapiens.GRCh38.91.gtf \
  --no_annotation -resume \
  -with-report report_s2.html
```

Because the annotation setting feeds only the `BAMBU` process, `-resume` caches
alignment, sorting and QC, and only transcript discovery runs again.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--reads` | `fastq/*.fastq.gz` | Glob matching the input fastq files. Sample names come from the filenames. |
| `--genome` | `ref/...dna_sm.primary_assembly.fa` | Genome fasta. |
| `--gtf` | `ref/Homo_sapiens.GRCh38.91.gtf` | Reference annotation. |
| `--outdir` | `results` | Output directory. |
| `--no_annotation` | `false` | Run Bambu without annotations (de novo discovery). |
| `--ndr` | `0.5` | Novel discovery rate, used only when no annotation is supplied. |
| `--min_reads` | `100000` | QC threshold: minimum reads per sample. |
| `--min_mapped_pct` | `80` | QC threshold: minimum percentage of reads mapped. |

## Output layout

```
results/
├── bam/      sorted BAM files and .bai indexes
├── qc/       per-sample metrics and qc_summary.tsv
└── bambu/    extended_annotations.gtf, count tables, se.rds
```

## Repository contents

- `main.nf` — pipeline definition
- `nextflow.config` — conda environments and per-process resources
- `bin/run_bambu.R` — the R script called by the `BAMBU` process
- `task4_bambu.R` — standalone Bambu run for Task 4
- `task4_answers.R` — computes the Task 4 answers from the saved result

## References

- Li, H. (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34(18):3094–3100.
- Chen, Y. et al. (2023) Context-aware transcript quantification from long-read RNA-seq data with bambu. *Nature Methods* 20:1187–1195.
- Di Tommaso, P. et al. (2017) Nextflow enables reproducible computational workflows. *Nature Biotechnology* 35:316–319.
- Danecek, P. et al. (2021) Twelve years of SAMtools and BCFtools. *GigaScience* 10(2).
- Chen, Y. et al. (2025) A systematic benchmark of Nanopore long-read RNA sequencing for transcript-level analysis in human cell lines. *Nature Methods*.
