test_that("Long test of Read_SBS_VCF_and_BAMs_to_verify_DBSs on input/HepG2_AA1.vcf -- full set of 636 DBSs", {
  skip_if(Sys.getenv("DO_LONG_TEST") != "y")
  Sys.time()
  stime <- system.time(
    xx <- Read_SBS_VCF_and_BAMs_to_verify_DBSs(
      input.vcf   = "input/HepG2_AA1.vcf",
      Nbam.name        = "input/HepG2_AA1_DBSlocs_Normal.bam",
      Tbam.name        = "input/HepG2_AA1_DBSlocs_Tumor.bam",
      variant.caller   = "strelka",
      verbose          = 1,
      unlink.slice.dir = FALSE)
  )
  Sys.time()
  # cat("\nLong test timing:\n", names(stime), "\n", stime, "\n")
  # cat(round(sum(stime[1:2]) / nrow(xx$evaluated.vcf), digits = 3),
  #     "CPU seconds per DBS\n")

  new <- data.table::fread(xx$evaluated.vcf.name)
  old <- data.table::fread(paste0(xx$evaluated.vcf.name, ".regress"))
  unlink(xx$evaluated.vcf.name)
  expect_equal(old, new)
  RegressSAMDirectory(old.dir = "input/long.test.N.slice.dir.regress/",
                      new.dir = xx$N.slice.dir)
  RegressSAMDirectory(old.dir = "input/long.test.T.slice.dir.regress/",
                      new.dir = xx$T.slice.dir)
  unlink(xx$N.slice.dir, recursive = TRUE)
  unlink(xx$T.slice.dir, recursive = TRUE)
})
