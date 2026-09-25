suppressPackageStartupMessages({library(bambu); library(SummarizedExperiment); library(GenomicRanges)})
se <- readRDS("bambu_out/se.rds")
rd <- rowData(se); cts <- assays(se)$counts
sink("bambu_out/task4_answers.txt", split = TRUE)
cat("Assays:", paste(assayNames(se), collapse=", "), "\n")
cat("rowData:", paste(colnames(rd), collapse=", "), "\n")
gtf <- rtracklayer::import("ref/Homo_sapiens.GRCh38.91.gtf")
gn <- unique(as.data.frame(mcols(gtf)[gtf$type=="gene", c("gene_id","gene_name")]))
sym <- function(i) gn$gene_name[match(i, gn$gene_id)]
rm(gtf); invisible(gc())
novel <- if ("novelTranscript" %in% colnames(rd)) rd$novelTranscript else grepl("^BambuTx", rownames(se))
cat("\n== 4.2a novel transcripts:", sum(novel), "\n")
sg <- transcriptToGeneExpression(se); gc5 <- sort(rowSums(assays(sg)$counts), decreasing=TRUE)[1:5]
cat("\n== 4.2b top 5 genes ==\n"); print(data.frame(gene=names(gc5), name=sym(names(gc5)), counts=round(gc5)), row.names=FALSE)
tpg <- table(rd$GENEID)
cat("\n== 4.2c all: min", min(tpg), "max", max(tpg), "mean", round(mean(tpg),3), "\n")
keep <- rowSums(cts) >= 10; t2 <- table(rd$GENEID[keep])
cat("== 4.2d expressed(>=10 total): transcripts", sum(keep), "genes", length(t2),
    "| min", min(t2), "max", max(t2), "mean", round(mean(t2),3), "\n")
ne <- lengths(rowRanges(se))
cat("== 4.2e single-exon novel:", sum(novel & ne==1), "\n")
if ("uniqueCounts" %in% assayNames(se)) {
  ou <- rowSums(cts) > 0 & rowSums(abs(cts - assays(se)$uniqueCounts)) < 1e-6
  tot <- sort(rowSums(cts)[ou], decreasing=TRUE)
  cat("== 4.2f unique-only transcripts:", sum(ou), "| top:\n"); print(round(head(tot,5)))
}
cat("\n== 4.4 candidates ==\n")
ids <- rownames(se)[novel]; rg <- unlist(range(rowRanges(se)[ids]))
cand <- data.frame(tx=ids, gene_name=sym(rd$GENEID[novel]), exons=ne[novel],
                   counts=round(rowSums(cts)[novel],1), samples=rowSums(cts[novel,,drop=FALSE]>=1),
                   ucsc=paste0("chr", seqnames(rg), ":", start(rg), "-", end(rg)))
if ("txClassDescription" %in% colnames(rd)) cand$class <- rd$txClassDescription[novel]
print(head(cand[order(-cand$counts),], 20), row.names=FALSE)
cat("\nBambu version:", as.character(packageVersion("bambu")), "\n")
sink()
