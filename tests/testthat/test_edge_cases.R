test_that("No read support at all", {
  xx <- DBSverify::ReadVCFAndBAMsAndVerifyDBSs(
    input.vcf.name = "input/1-30.vcf",
    Nbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
    Tbam.name = "input/HepG2_AA1_DBSlocs_Normal.bam",
    variant.caller = "strelka",
    num.cores = 10,
    N.slice.dir = "input/1-30_test_N",
    T.slice.dir = "input/1-30_test_T",
    unlink.slice.dir = TRUE,
    verbose = 0)
  new <- data.table::fread(xx$evaluated.vcf.name)
  old <- data.table::fread(paste0(xx$evaluated.vcf.name, ".regress"))
  unlink(xx$evaluated.vcf.name)
  expect_equal(old, new)
})

test_that("No DBSs at all", {
  xx <- DBSverify::ReadVCFAndBAMsAndVerifyDBSs(
    input.vcf.name = "input/no-dbs.vcf",
    Nbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
    Tbam.name = "input/HepG2_AA1_DBSlocs_Normal.bam",
    variant.caller = "strelka",
    num.cores = 10,
    N.slice.dir = "input/1-30_test_N",
    T.slice.dir = "input/1-30_test_T",
    unlink.slice.dir = TRUE,
    verbose = 0)
  new <- data.table::fread(xx$evaluated.vcf.name)
  old <- data.table::fread(paste0(xx$evaluated.vcf.name, ".regress"))
  unlink(xx$evaluated.vcf.name)
  expect_equal(old, new)
})

