test_that(
  "No normal BAM index",
  expect_error(
    xx <- Read_SBS_VCF_and_BAMs_to_verify_DBSs(
      input.vcf = "input/6-54496940.vcf",
      Nbam.name = "input/bam_no_bai.bam",
      Tbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
      variant.caller = "strelka"),
    regexp = "BAM index file"
  ))

test_that(
  "No tumor BAM index",
  expect_error(
    xx <- Read_SBS_VCF_and_BAMs_to_verify_DBSs(
      input.vcf = "input/6-54496940.vcf",
      Nbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
      Tbam.name = "input/bam_no_bai.bam",
      variant.caller = "strelka"),
    regexp = "BAM index file"
  ))

test_that(
  "No VCF file",
  expect_error(
    xx <- Read_SBS_VCF_and_BAMs_to_verify_DBSs(
      input.vcf = "nope",
      Nbam.name = "input/HepG2_AA1_DBSlocs_Normal.bam",
      Tbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
      variant.caller = "strelka"),
    regexp = "does not exist"
  ))

test_that(
  "No normal BAM file",
  expect_error(
    xx <- Read_SBS_VCF_and_BAMs_to_verify_DBSs(
      input.vcf = "input/6-54496940.vcf",
      Nbam.name = "nope",
      Tbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
      variant.caller = "strelka"),
    regexp = "BAM file"
  ))

test_that(
  "No tumor BAM file",
  expect_error(
    xx <- Read_SBS_VCF_and_BAMs_to_verify_DBSs(
      input.vcf = "input/6-54496940.vcf",
      Nbam.name = "input/HepG2_AA1_DBSlocs_Normal.bam",
      Tbam.name = "nope",
      variant.caller = "strelka"),
    regexp = "BAM file"
  ))

test_that(
  "N.slice.dir == T.slice.dir",
  expect_error(
    xx <- Read_SBS_VCF_and_BAMs_to_verify_DBSs(
      input.vcf = "input/6-54496940.vcf",
      Nbam.name = "input/HepG2_AA1_DBSlocs_Normal.bam",
      Tbam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam",
      variant.caller = "strelka",
      T.slice.dir = "tmp.dir",
      N.slice.dir = "tmp.dir"),
    regexp = "must be different"
  ))
