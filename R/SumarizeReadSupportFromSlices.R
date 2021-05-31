SummarizeReadSupportFromSlices <- function(vcf, N.slice.dir, T.slice.dir) {
  for(v in 1:nrow(vcf)){

    CHROM <- vcf[v, "CHROM"]
    POS   <-vcf[v, "POS"]
    REF   <- vcf[v, "REF"]
    ALT   <- vcf[v, "ALT"]

    vcf[v, "NreadSupport"] <-
      Slice2ReadSupport(N.slice.dir, CHROM = CHROM, POS = POS, REF = REF, ALT = ALT)

    vcf[v, "TreadSupport"] <-
      Slice2ReadSupport(T.slice.dir, CHROM = CHROM, POS = POS, REF = REF, ALT = ALT)
  }
  return(vcf)
}
