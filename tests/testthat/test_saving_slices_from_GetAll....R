test_that("Test ReadVCFAndBAMsAndVerifyDBSs BAM slices in tempdir()", {
  xx <- ReadVCFAndBAMsAndVerifyDBSs(
      vcf.name         = "input/nine.vcf",
      Nbam.name        = "input/HepG2_AA1_DBSlocs_Normal.bam",
      Tbam.name        = "input/HepG2_AA1_DBSlocs_Tumor.bam",
      variant.caller   = "strelka",
      num.cores        = 1,
      unlink.slice.dir = FALSE)
    new <- data.table::fread(xx$vcf.name)
    old <- data.table::fread(paste0(xx$vcf.name, ".regress"))
    unlink(xx$vcf.name)
    expect_equal(old, new)
    RegressSAMDirectory(old.dir = "input/nine.regress.normal", new.dir = xx$N.slice.dir)
    RegressSAMDirectory(old.dir = "input/nine.regress", new.dir = xx$T.slice.dir)
})

test_that("Test ReadVCFAndBAMsAndVerifyDBSs saving BAM slices in sepcified folder", {
  xx <- ReadVCFAndBAMsAndVerifyDBSs(
    vcf.name         = "input/nine.vcf",
    Nbam.name        = "input/HepG2_AA1_DBSlocs_Normal.bam",
    Tbam.name        = "input/HepG2_AA1_DBSlocs_Tumor.bam",
    variant.caller   = "strelka",
    num.cores        = 1,
    N.slice.dir      = "tmp.test.N.slice.dir",
    T.slice.dir      = "tmp.test.T.slice.dir",
    unlink.slice.dir = FALSE)
  new <- data.table::fread(xx$vcf.name)
  old <- data.table::fread(paste0(xx$vcf.name, ".regress"))
  unlink(xx$vcf.name)
  expect_equal(old, new)
  RegressSAMDirectory(old.dir = "input/nine.regress.normal", new.dir = xx$N.slice.dir)
  RegressSAMDirectory(old.dir = "input/nine.regress", new.dir = xx$T.slice.dir)
  unlink("tmp.test.N.slice.dir", recursive = TRUE)
  unlink("tmp.test.T.slice.dir", recursive = TRUE)
})
