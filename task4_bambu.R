## DSA4262 Assignment 1 - Task 4: transcript discovery & quantification with Bambu
suppressPackageStartupMessages({
  library(bambu)
  library(SummarizedExperiment)
  library(GenomicRanges)
})

fa.file   <- "ref/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
gtf.file  <- "ref/Homo_sapiens.GRCh38.91.gtf"
bam.dir   <- "bam"
n.cores   <- 4
min.reads <- 10

bam.files <- list.files(bam.dir, pattern = "\\.bam$", full.names = TRUE)
stopifnot(length(bam.files) > 0, file.exists(fa.file), file.exists(gtf.file))
dir.create("bambu_out", showWarnings = FALSE)

if (!file.exists("bambu_out/se.rds")) {
  annotations <- prepareAnnotations(gtf.file)
  se <- bambu(reads = bam.files, annotations = annotations,
              genome = fa.file, ncore = n.cores)
  saveRDS(se, "bambu_out/se.rds")
  writeBambuOutput(se, path = "bambu_out/")
} else {
  se <- readRDS("bambu_out/se.rds")
}

sink("bambu_out/task4_answers.txt", split = TRUE)

cat("\n===== 4.1 show(se) =====\n")
show(se)

rd  <- rowData(se)
cts <- assays(se)$counts
cat("\nrowData columns:", paste(colnames(rd), collapse = ", "), "\n")
cat("Assays:", paste(assayNames(se), collapse = ", "), "\n")

gtf <- rtracklayer::import(gtf.file)
gene.names <- unique(as.data.frame(mcols(gtf)[gtf$type == "gene",
                                              c("gene_id", "gene_name")]))
sym <- function(ids) gene.names$gene_name[match(ids, gene.names$gene_id)]
rm(gtf); invisible(gc())

novel <- if ("novelTranscript" %in% colnames(rd)) rd$novelTranscript else
  grepl("^BambuTx", rownames(se))
cat("\n===== 4.2a =====\n")
cat("Novel transcripts:", sum(novel), "\n")
cat("(Novel genes:", length(unique(rd$GENEID[grepl("^BambuGene", rd$GENEID)])), ")\n")

se.gene <- transcriptToGeneExpression(se)
g.counts <- rowSums(assays(se.gene)$counts)
top5 <- head(sort(g.counts, decreasing = TRUE), 5)
cat("\n===== 4.2b (total counts across all samples) =====\n")
print(data.frame(gene_id = names(top5), gene_name = sym(names(top5)),
                 total_counts = round(top5)), row.names = FALSE)
if ("CPM" %in% assayNames(se.gene)) {
  g.cpm <- rowMeans(assays(se.gene)$CPM)
  top5c <- head(sort(g.cpm, decreasing = TRUE), 5)
  cat("\n(Alternative: mean CPM, adjusts for sequencing depth)\n")
  print(data.frame(gene_id = names(top5c), gene_name = sym(names(top5c)),
                   mean_CPM = round(top5c, 1)), row.names = FALSE)
}

tx.per.gene <- table(rd$GENEID)
cat("\n===== 4.2c (extended annotation, all transcripts) =====\n")
cat("Min:", min(tx.per.gene), " Max:", max(tx.per.gene),
    " Mean:", round(mean(tx.per.gene), 3), "\n")
cat("Gene with max:", names(which.max(tx.per.gene)),
    sym(names(which.max(tx.per.gene))), "\n")

expr.total <- rowSums(cts) >= min.reads
expr.any   <- rowSums(cts >= min.reads) >= 1
for (nm in c("total", "any")) {
  keep <- if (nm == "total") expr.total else expr.any
  t2 <- table(rd$GENEID[keep])
  cat("\n===== 4.2d (", min.reads, "+ reads,",
      ifelse(nm == "total", "summed across samples", "in at least one sample"),
      ") =====\n")
  cat("Expressed transcripts:", sum(keep), " Genes:", length(t2), "\n")
  cat("Min:", min(t2), " Max:", max(t2), " Mean:", round(mean(t2), 3), "\n")
}

n.exons <- lengths(rowRanges(se))
cat("\n===== 4.2e =====\n")
cat("Single-exon novel transcripts:", sum(novel & n.exons == 1), "\n")

cat("\n===== 4.2f =====\n")
if ("uniqueCounts" %in% assayNames(se)) {
  uni <- assays(se)$uniqueCounts
  only.unique <- rowSums(cts) > 0 & rowSums(abs(cts - uni)) < 1e-6
  tot <- rowSums(cts)[only.unique]
  best <- names(which.max(tot))
  cat("Transcripts whose counts are all unique reads:", sum(only.unique), "\n")
  cat("Highest:", best, "| gene", rd[best, "GENEID"], sym(rd[best, "GENEID"]),
      "| total counts", round(max(tot)), "\n")
  print(head(sort(tot, decreasing = TRUE), 5))
} else {
  cat("No 'uniqueCounts' assay - check assayNames(se) above\n")
}

cat("\n===== 4.4 candidate novel transcripts =====\n")
ids <- rownames(se)[novel]
rg  <- unlist(range(rowRanges(se)[ids]))
chr <- as.character(seqnames(rg)); chr[chr == "MT"] <- "M"
cand <- data.frame(tx = ids, gene = rd$GENEID[novel], gene_name = sym(rd$GENEID[novel]),
                   exons = n.exons[novel], strand = as.character(strand(rg)),
                   total_counts = round(rowSums(cts)[novel], 1),
                   samples_with_reads = rowSums(cts[novel, , drop = FALSE] >= 1),
                   ucsc = paste0("chr", chr, ":", start(rg), "-", end(rg)))
if ("txClassDescription" %in% colnames(rd)) cand$class <- rd$txClassDescription[novel]
print(head(cand[order(-cand$total_counts), ], 20), row.names = FALSE)

cat("\n===== 4.3 versions =====\n")
cat("Bambu version:", as.character(packageVersion("bambu")), "\n")
cat("Samples:", paste(colnames(se), collapse = ", "), "\n")
cat("Metadata fields:", paste(names(metadata(se)), collapse = ", "), "\n")
print(sessionInfo())

sink()
cat("\nDone - answers saved to bambu_out/task4_answers.txt\n")
