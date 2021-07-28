SummarizeReadSupportFromSlices <- function(vcf, N.slice.dir, T.slice.dir) {
  for(v in 1:nrow(vcf)){

    CHROM <- vcf[v, "CHROM"]
    POS   <-vcf[v, "POS"]
    REF   <- vcf[v, "REF"]
    ALT   <- vcf[v, "ALT"]

    rsN <- Slice2ReadSupport(N.slice.dir,
                             CHROM = CHROM, POS = POS, REF = REF, ALT = ALT)
    vcf[v, "NreadSupport"] <- rsN$read.support

    rsT <- Slice2ReadSupport(T.slice.dir,
                             CHROM = CHROM, POS = POS, REF = REF, ALT = ALT)
    vcf[v, "TreadSupport"]             <- rsT$read.support
    vcf[v, "num_bad_mapped_reads"]     <- rsT$num.bad.mapped.reads
    vcf[v, "num_bad_mapped_DBS_reads"] <- rsT$num.bad.mapped.DBS.reads

  }
  return(vcf)
}
