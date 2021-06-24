Remove_non_canonical_chromosomes <- function(vcf, verbose = 0) {
  stopifnot(grepl("(#|)CHROM$", colnames(vcf)[1]))
  ok.rows <- grep("^(chr|)((\\d+)|X|x|Y|y)$", vcf[ , 1])
  if (verbose && (length(ok.rows) < nrow(vcf))) {
    bad.rows <- vcf[-ok.rows, ]
    bad.chr <- unique(bad.rows[ , 1])
    message("Removing ", nrow(bad.rows),
            " rows with the following noncanonical chromosomes: ",
            paste(bad.chr, collapse = ", "))
  }
  return(vcf[ok.rows, ])
}
