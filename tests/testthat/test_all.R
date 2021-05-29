test_that("Long test of input/HepG2_AA1_20uM_SL_cl1_SNVresult.vcf -- full set of DBSs", {
skip_if(Sys.getenv("DO_LONG_TEST") == "")
xx <- ReadVCFAndBAMsAndProcess(
  vcf.name  = "input/HepG2_AA1_20uM_SL_cl1_SNVresult.vcf",
  Nbam.name = "input/HepG2_AA1_DBSlocs_Normal.bam",
  Tbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
  variant.caller = "strelka",
  num.cores      = 1)

new <- data.table::fread(xx$vcf.name)
old <- data.table::fread(paste0(xx$vcf.name, ".regress"))
unlink(xx$vcf.name)
expect_equal(old, new)

})
