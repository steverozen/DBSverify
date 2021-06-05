basic.test <- function(input.vcf.root) {
  vcf.name <- paste0("input/", input.vcf.root, ".vcf")
  xx <- ReadVCFAndBAMsAndVerifyDBSs(
    input.vcf.name   = vcf.name,
    Nbam.name        = "input/HepG2_AA1_DBSlocs_Normal.bam",
    Tbam.name        = "input/HepG2_AA1_DBSlocs_Tumor.bam",
    variant.caller   = "strelka",
    num.cores        = 1,
    unlink.slice.dir = TRUE)
  new <- data.table::fread(xx$evaluated.vcf.name)
  old <- data.table::fread(paste0(xx$evaluated.vcf.name, ".regress"))
  unlink(xx$evaluated.vcf.name)
  expect_equal(old, new)
}

# example bash code to generate input VCF
# cd tests/testthat/input
# (cat vcf.header; grep -A1 54496939 SNVresult.vcf) > 6-54496940.vcf
#  Originally HepG2_AA1_20uM_SL_cl1_SNVresult.vcf

test_that("True DBS",
          basic.test("1-74823446"))

test_that("Adjacent SBSs - 1st pos only in some reads, DBS in other reads",
          basic.test("1-116218637"))

test_that("Adjacent SBSs - opposite alleles",
          basic.test("5-100759123"))

test_that("1st position overlaps germline SNP",
          basic.test("6-54496940"))

test_that("Germline DBS",
          basic.test("6-98489243"))

test_that("Adjacent SBSs",
          basic.test("7-131714192")) # Tumor reads supporting DBS and only supporting position 2 in germline and tumor

test_that("Adjacent SBSs",
          basic.test("8-106268966")) # This actually only has reads supporting 1 position

test_that("Adjacent SBSs - 2nd pos only in some reads, DBS in other reads",
          basic.test("10-59377023"))

test_that("2nd position overlaps germline SNP",
          basic.test("15-94501219"))



