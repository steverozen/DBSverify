if (FALSE) {
  foo <- HMF_prep_BEDs_and_DBS_VCFs(input.vcf.dir = "~/mvv/HMF-VCFs/")

  # This was the command used to generate the DBS-only VCF
  # and BED files for Hartwig Medication Foundation (HMF)
  # data.
  foo <- HMF_prep_BEDs_and_DBS_VCFs("~/mvv/HMF_all_input_VCF/",
                                    out.bed.dir = "~/mvv/HMF_all_BED_and_VCF",
                                    out.vcf.dir = "~/mvv/HMF_all_BED_and_VCF")
}

HMF_prep_BEDs_and_DBS_VCFs <- function(input.vcf.dir,
                                       out.bed.dir = input.vcf.dir,
                                       out.vcf.dir = input.vcf.dir,
                                       verbose = TRUE) {

  if (verbose) message("HMV VCFs from ", input.vcf.dir)

  input.vcfs <- dir(input.vcf.dir,
                    full.names = FALSE,
                    # pattern = "CPCT02030256T.purple.somatic.vcf.gz")
                    pattern = ".vcf.gz")

  # ignore possible index files

  input.vcfs <-
    grep("(\\.idx$)|(\\.tbi$)", input.vcfs, value = TRUE, invert = TRUE)

  xx <- lapply(input.vcfs,
               HMF_prep_1_BED_and_DBS_VCF,
               out.bed.dir = out.bed.dir,
               out.vcf.dir = out.vcf.dir,
               verbose     = verbose,
               input.vcf.dir = input.vcf.dir)
  summary <- data.frame(t(matrix(unlist(xx), nrow = 3)))
  return(invisible(summary))
}
