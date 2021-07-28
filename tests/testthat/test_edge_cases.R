test_that("No read support at all from SBS VCF 1", {
  xx <- DBSverify::Read_SBS_VCF_and_BAMs_to_verify_DBSs(
    input.vcf = "input/1-30.vcf",
    Nbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
    Tbam.name = "input/HepG2_AA1_DBSlocs_Normal.bam",
    variant.caller = "strelka",
    N.slice.dir = "input/1-30_test_N",
    T.slice.dir = "input/1-30_test_T",
    unlink.slice.dir = TRUE,
    verbose = 0)
  new <- data.table::fread(xx$evaluated.vcf.name)
  old <- data.table::fread(paste0(xx$evaluated.vcf.name, ".regress"))
  unlink(xx$evaluated.vcf.name)
  expect_equal(old, new)
})

test_that("No read support at all from DBS VCF 2", {
  xx <- DBSverify::Read_DBS_VCF_and_BAMs_to_verify_DBSs(
    input.vcf = "input/1-30-DBS.vcf",
    Nbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
    Tbam.name = "input/HepG2_AA1_DBSlocs_Normal.bam",
    N.slice.dir = "input/1-30_test_N",
    T.slice.dir = "input/1-30_test_T",
    unlink.slice.dir = TRUE,
    verbose = 0)
  new <- data.table::fread(xx$evaluated.vcf.name)
  old <- data.table::fread(paste0(xx$evaluated.vcf.name, ".regress"))
  unlink(xx$evaluated.vcf.name)
  expect_equal(old, new)
})

test_that("No DBSs at all from SBS VCF", {
  xx <- DBSverify::Read_SBS_VCF_and_BAMs_to_verify_DBSs(
    input.vcf = "input/no-dbs.vcf",
    Nbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
    Tbam.name = "input/HepG2_AA1_DBSlocs_Normal.bam",
    variant.caller = "strelka",
    N.slice.dir = "input/1-30_test_N",
    T.slice.dir = "input/1-30_test_T",
    unlink.slice.dir = TRUE,
    verbose = 0)
  new <- data.table::fread(xx$evaluated.vcf.name)
  old <- data.table::fread(paste0(xx$evaluated.vcf.name, ".regress"))
  unlink(xx$evaluated.vcf.name)
  expect_equal(old, new)
})

