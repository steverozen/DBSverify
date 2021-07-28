test_that("Test GetAllBAMSlicesSamtools use tempfile", {
  aa <- ICAMS::ReadAndSplitVCFs(# "input/10-59377023.vcf",
    "input/nine.vcf",
    # Nine examples of various different DBSconclusion and read support results
    variant.caller = "strelka",
    max.vaf.diff   = 1,
    num.cores      = 1)
  dbs.vcf <- aa$DBS[[1]]
  rm(aa)

  tmp.output <- "tmp.slice.output"
  for (spec.outdir in c(TRUE, FALSE)) {
    if (spec.outdir) {
      stime <- system.time(
        outdir <- GetAllBAMSlicesSamtools(
          where.to.put.slices = tmp.output,
          vcf = dbs.vcf,
          bam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam"))
    } else {
      stime <- system.time(
        outdir <- GetAllBAMSlicesSamtools(
          vcf = dbs.vcf,
          bam.name = "input/HepG2_AA1_DBSlocs_Tumor.bam"))
    }

    # cat("\nTest timing:\n", names(stime), "\n", stime, "\n")
    # cat(round(sum(stime[1:2]) / nrow(dbs.vcf), digits = 3),
    #    "CPU seconds per DBS\n")

    sams <- dir(outdir, pattern = ".sam")
    for (ss in sams) {
      new <- ReadSamfile(file.path(outdir, ss))
      old <- ReadSamfile(file.path("input/nine.regress/", ss))
      expect_equal(old, new)
    }
  } # for

  unlink(tmp.output, recursive = TRUE)

})
