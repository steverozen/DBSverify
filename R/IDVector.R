IDVector <- function(tmp.vcf) {
  apply(tmp.vcf[ , c("CHROM", "POS", "REF", "ALT"), drop =  FALSE],
        MARGIN = 1,
        function(x) {
          rr <- paste(x, collapse = ":")
          return(gsub(" ", "", rr, fixed = TRUE))
        })
}
