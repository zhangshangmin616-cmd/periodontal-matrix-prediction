# Recreate barcode-level scDblFinder calls for one QC-filtered sparse count matrix.
# Usage: Rscript run_doublets.R <library> <cache-directory> <output-directory>
suppressPackageStartupMessages({
  library(Matrix)
  library(scDblFinder)
})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
key <- args[1L]
cache_dir <- args[2L]
dest <- args[3L]
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
part <- function(suffix) file.path(cache_dir, paste0(key, "_dbl_", suffix))
read_int32 <- function(path) {
  bytes <- file.size(path)
  stopifnot(bytes %% 4 == 0)
  readBin(path, integer(), n = bytes / 4, size = 4L, endian = "little")
}
genes <- readLines(part("genes.tsv"), warn = FALSE)
barcodes <- readLines(part("barcodes.tsv"), warn = FALSE)
stopifnot(!anyDuplicated(genes), !anyDuplicated(barcodes), all(nzchar(barcodes)))
m <- new("dgCMatrix", i = read_int32(part("i.bin")), p = read_int32(part("p.bin")),
         x = as.numeric(read_int32(part("x.bin"))),
         Dim = as.integer(c(length(genes), length(barcodes))), Dimnames = list(genes, barcodes))
stopifnot(validObject(m), min(m@x) > 0)
set.seed(20261124)
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
rm(m)
invisible(gc())
print(Sys.time())
sce <- scDblFinder(sce, clusters = FALSE, dbr.per1k = 0.008,
                   BPPARAM = BiocParallel::SerialParam(RNGseed = 20261124, progressbar = FALSE))
res <- data.frame(barcode = colnames(sce), score = sce$scDblFinder.score,
                  call = as.character(sce$scDblFinder.class))
stopifnot(identical(res$barcode, barcodes), nrow(res) == length(barcodes),
          all(is.finite(res$score)), all(res$call %in% c("doublet","singlet")))
tmp <- file.path(dest, paste0(key,".tsv.tmp"))
write.table(res, file = tmp, quote = FALSE, row.names = FALSE, sep = "\t")
path <- file.path(dest, paste0(key,".tsv"))
if (!file.rename(tmp, path)) stop("Could not write barcode results")
print(table(res$call))
print(Sys.time())
writeLines(capture.output(sessionInfo()), file.path(dest, paste0(key, "_session.txt")))
