debug(VerifyDBSVcf)
debug(GetAllBAMSlicesSamtools)

xx <- DBSverify::ReadVCFAndBAMsAndVerifyDBSs(
  # input.vcf.name = "../mvv/CPCT02020306T.purple.somatic.vcf.gz",
  input.vcf.name = "../mvv/WIDE01010310T.purple.somatic.vcf.gz", # this one is short
  # input.vcf.name = "../mvv/CPCT02030256T.purple.somatic.vcf.gz",
  Nbam.name = "tests/testthat/input/HepG2_AA1_DBSlocs_Tumor.bam",
  Tbam.name = "tests/testthat/input/HepG2_AA1_DBSlocs_Normal.bam",
  variant.caller = "unknown",
  num.cores = 10,
  N.slice.dir = "../mvv/HMF_test_N",
  T.slice.dir = "../mvv/HMF_test_T",
  unlink.slice.dir = FALSE,
  verbose = 10
)

