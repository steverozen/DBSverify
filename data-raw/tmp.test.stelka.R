foo <-
  Read_SBS_VCF_and_BAMs_to_verify_DBSs(
    input.vcf = "~/mvv/UTUC-test-case/passed.somatic.snvs.vcf.gz",
    Nbam.name = "~/mvv/UTUC-test-case/corrected_AA_UTUC_Taiwan_N44.rmdup.bam",
    Tbam.name = "~/mvv/UTUC-test-case/corrected_AA_UTUC_Taiwan_T11.rmdup.bam",
    variant.caller = "strelka",
    outfile = "~/mvv/tmp.test.strelka/tmp.test.strelka2.vcf",
    verbose = 1

  )

v1 <- data.table::fread("~/mvv/tmp.test.strelka/tmp.test.strelka.vcf")
v2 <- data.table::fread("~/mvv/tmp.test.strelka/tmp.test.strelka2.vcf")
colnames(v2)[3:18] <- paste(colnames(v2)[3:18], "_v2")
vx <- dplyr::full_join(v1, v2)
