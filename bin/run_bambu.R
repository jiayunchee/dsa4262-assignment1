#!/usr/bin/env Rscript
## Called by the BAMBU process in main.nf
## Usage: run_bambu.R <genome.fa> <annotation.gtf|NONE> <ndr> <ncore> <bam1> <bam2> ...

args   <- commandArgs(trailingOnly = TRUE)
genome <- args[1]
gtf    <- args[2]
ndr    <- as.numeric(args[3])
ncore  <- as.integer(args[4])
bams   <- args[-(1:4)]

suppressPackageStartupMessages(library(bambu))

if (identical(gtf, "NONE")) {
  ## Scenario 2: no annotations. Bambu cannot estimate an NDR threshold from an
  ## annotation here, so one is supplied explicitly.
  se <- bambu(reads = bams, annotations = NULL, genome = genome,
              ncore = ncore, NDR = ndr)
} else {
  ## Scenario 1: discovery and quantification against the supplied annotation
  annotations <- prepareAnnotations(gtf)
  se <- bambu(reads = bams, annotations = annotations, genome = genome,
              ncore = ncore)
}

saveRDS(se, "se.rds")
writeBambuOutput(se, path = ".")
show(se)
