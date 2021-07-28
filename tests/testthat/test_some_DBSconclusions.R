jig <- function(old, new) {
  View(dplyr::full_join(old, new))
  View(new)
}

basic.test <- function(input.vcf.root, unlink.me = TRUE) {
  vcf.name <- paste0("input/", input.vcf.root, ".vcf")
  xx <- Read_SBS_VCF_and_BAMs_to_verify_DBSs(
    input.vcf        = vcf.name,
    Nbam.name        = "input/HepG2_AA1_DBSlocs_Normal.bam",
    Tbam.name        = "input/HepG2_AA1_DBSlocs_Tumor.bam",
    variant.caller   = "strelka",
    unlink.slice.dir = TRUE)
  new <- data.table::fread(xx$evaluated.vcf.name)
  old <- data.table::fread(paste0(xx$evaluated.vcf.name, ".regress"))
  if (unlink.me) unlink(xx$evaluated.vcf.name)
  # jig(old, new)
  expect_equal(old, new)
}

# example bash code to generate input VCF
# cd tests/testthat/input
# (cat vcf.header; grep -A1 54496939 SNVresult.vcf) > 6-54496940.vcf
#  Originally HepG2_AA1_20uM_SL_cl1_SNVresult.vcf

test_that("True DBS",
          basic.test("1-74823446")) # ok

test_that("True DBS, requires discarding some reads",
          basic.test("1-116218637")) # ok

test_that("Adjacent SBSs - opposite alleles",
          basic.test("5-100759123")) # ok

test_that("Neither position supported",
          basic.test("6-54496940")) # ok

test_that("True DBS",
          basic.test("6-98489243")) # ok

test_that("True DBS",
          basic.test("7-131714192")) # ok

test_that("Neither position supported",
          basic.test("8-106268966")) # ok Possibly investigate -- lots of reads discarded because of CIGAR

test_that("Adjacent SBSs - 2nd pos only in some reads, DBS in other reads",
          basic.test("10-59377023")) # ok

test_that("2nd position overlaps germline SNP",
          basic.test("15-94501219")) #xxxxxx bad
