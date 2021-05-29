test_that("Test 1 true DBS", {

xx <- ReadVCFAndBAMsAndProcess(
  vcf.name = "input/74823446.vcf",
  Nbam.name = "input/HepG2_AA1_DBSlocs_Normal.bam",
  Tbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
  variant.caller = "strelka",
  num.cores      = 1)

  new <- data.table::fread(xx$vcf.name)
  old <- data.table::fread(paste0(xx$vcf.name, ".regress"))
  unlink(xx$vcf.name)
  expect_equal(old, new)
})

