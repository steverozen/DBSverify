test_that(
  "No normal BAM index",
  expect_error(
    xx <- ReadVCFAndBAMsAndVerifyDBSs(
      vcf.name = "input/6-54496940.vcf",
      Nbam.name = "input/bam_no_bai.bam",
      Tbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
      variant.caller = "strelka",
      num.cores      = 1),
    regexp = "BAM index file"
  ))

test_that(
  "No tumor BAM index",
  expect_error(
    xx <- ReadVCFAndBAMsAndVerifyDBSs(
      vcf.name = "input/6-54496940.vcf",
      Nbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
      Tbam.name = "input/bam_no_bai.bam",
      variant.caller = "strelka",
      num.cores      = 1),
    regexp = "BAM index file"
  ))

test_that(
  "No VCF file",
  expect_error(
    xx <- ReadVCFAndBAMsAndVerifyDBSs(
      vcf.name = "nope",
      Nbam.name = "input/HepG2_AA1_DBSlocs_Normal.bam",
      Tbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
      variant.caller = "strelka",
      num.cores      = 1),
    regexp = "does not exist"
  ))

test_that(
  "No normal BAM file",
  expect_error(
    xx <- ReadVCFAndBAMsAndVerifyDBSs(
      vcf.name = "input/6-54496940.vcf",
      Nbam.name = "nope",
      Tbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
      variant.caller = "strelka",
      num.cores      = 1),
    regexp = "BAM file"
  ))

test_that(
  "No tumor BAM file",
  expect_error(
    xx <- ReadVCFAndBAMsAndVerifyDBSs(
      vcf.name = "input/6-54496940.vcf",
      Nbam.name = "input/HepG2_AA1_DBSlocs_Normal.bam",
      Tbam.name = "nope",
      variant.caller = "strelka",
      num.cores      = 1),
    regexp = "BAM file"
  ))
