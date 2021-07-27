test_that("Is the IGV script ok", {
  ss <- test_create_IGV_batch_snapshot(run.IGV = FALSE,
                                       vcf.name = "input/HepG2_AA1_DBS_evaluated.vcf")
  xx <- readLines(con = "input/test_IGV_script.txt")
  expect_equal(xx[c(1:2,4:4462)], ss[c(1:2,4:4462)]) # Line 3 is a tempfile().
})
