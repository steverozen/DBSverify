test_that("Long test of input/HepG2_AA1.vcf -- full set of 636 DBSs", {
  skip_if(Sys.getenv("DO_LONG_TEST") == "")
  Sys.time()
  stime <- system.time(
    xx <- ReadVCFAndBAMsAndVerifyDBSs(
      input.vcf.name   = "input/HepG2_AA1.vcf",
      Nbam.name        = "input/HepG2_AA1_DBSlocs_Normal.bam",
      Tbam.name        = "input/HepG2_AA1_DBSlocs_Tumor.bam",
      variant.caller   = "strelka",
      num.cores        = 1,
      N.slice.dir      = "long.test.tmp.N.slice.dir",
      T.slice.dir      = "long.test.tmp.T.slice.dir",
      unlink.slice.dir = FALSE)
  )
  Sys.time()
  cat("\nLong test timing:\n", names(stime), "\n", stime, "\n")
  cat(round(sum(stime[1:2]) / nrow(xx$evaluated.vcf), digits = 3),
      "CPU seconds per DBS\n")

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
