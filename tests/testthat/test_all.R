test_that("Long test of input/HepG2_AA1_20uM_SL_cl1_SNVresult.vcf -- full set of 636 DBSs", {
  skip_if(Sys.getenv("DO_LONG_TEST") == "")
  Sys.time()
  stime <- system.time(
    xx <- ReadVCFAndBAMsAndProcess(
      vcf.name  = "input/SNVresult.vcf",
      Nbam.name = "input/HepG2_AA1_DBSlocs_Normal.bam",
      Tbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
      variant.caller = "strelka",
      num.cores      = 1))
  Sys.time()
  cat("\nLong test timing:\n", names(stime), "\n", stime, "\n")
  cat(round(sum(stime[1:2]) / nrow(xx$evaluated.vcf), digits = 3),
      "CPU seconds per DBS\n")

  new <- data.table::fread(xx$vcf.name)
  old <- data.table::fread(paste0(xx$vcf.name, ".regress"))
  unlink(xx$vcf.name)
  expect_equal(old, new)
})
