vcf.name <- "input/74823436-74823456.vcf"

vcf.list <-ICAMS::ReadAndSplitVCFs(vcf.name,
                                   variant.caller = "strelka",
                                   num.of.cores = 10,
                                   names.of.VCFs = smp,
                                   max.vaf.diff = 1)

DBS.vcf <- vcf.list$DBS[[1]]

evaluated.vcf <- VerifyDBSVcf(DBS.vcf,
                              Nbam = "input/HepG2_AA1_1:74823436-74823456_N.bam",
                              Tbam = "input/HepG2_AA1_1:74823436-74823456_T.bam")
