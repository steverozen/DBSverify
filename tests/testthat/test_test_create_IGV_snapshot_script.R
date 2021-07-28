test_that("Is the IGV script ok", {
  ss <- test_create_IGV_snapshot_script(
    run.IGV = FALSE,
    vcf.name = "input/HepG2_AA1_DBS_evaluated.vcf")
  xx <- readLines(con = "input/test_IGV_script.txt")
  # Some lines point to files that have locations that vary:
  diffs <- which(xx != ss)
  if (length(diffs) > 0) cat("\n\ndiffs ", paste(diffs, collapse = ", "))
  else cat("\n\nno diffs")
  lines.to.check <- c(1:2,4:6,8,10:4462)
  expect_equal(xx[lines.to.check], ss[lines.to.check])
})
