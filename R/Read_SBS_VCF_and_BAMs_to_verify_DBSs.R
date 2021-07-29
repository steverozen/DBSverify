#' Aggressively merge SBSs into DBSs then determine whether sequencing reads in fact support DBSs.
#'
#' @param variant.caller One of \code{"strelka"}, \code{"PCAWG"},
#'  or \code{"unknown"}.
#'  Merging adjacent SBS is done by \code{ReadAndSplitVCFs} in the ICAMS package.
#'
#' @inheritParams Read_DBS_VCF_and_BAMs_to_verify_DBSs
#'
#' @details Note: argument \code{input.vcf} must be a file path.
#'  This function creates a new VCF file.
#'  See \code{\link{Read_DBS_VCF_and_BAMs_to_verify_DBSs}}
#'  for details.
#'
#' @return  Same as \code{\link{Read_DBS_VCF_and_BAMs_to_verify_DBSs}}.
#'
#' @md
#' @export
#'

Read_SBS_VCF_and_BAMs_to_verify_DBSs <- function(input.vcf,
                                                 Nbam.name,
                                                 Tbam.name,
                                                 variant.caller,
                                                 N.slice.dir      = tempfile(),
                                                 T.slice.dir      = tempfile(),
                                                 unlink.slice.dir = TRUE,
                                                 outfile          = NULL,
                                                 verbose          = 0) {

  CheckBAM(Nbam.name)
  CheckBAM(Tbam.name)
  TestBAMAndSamtools(Nbam.name)
  TestBAMAndSamtools(Tbam.name)

  if (!is.null(outfile)) {
    cat("test", file = outfile)
    unlink(outfile)
  }

  if (N.slice.dir == T.slice.dir) {
    stop("T.slice.dir and N.slice.dir must be different; got ",
         N.slice.dir)
  }

  if (variant.caller %in% c("strelka", "PCAWG", "unknown")) {
    get.vaf.function <- NULL
    if (variant.caller == "PCAWG") {
      variant.caller <- "unknown"
      # get.vaf.function <- GetPCAWGVAF
    }
  } else {
    stop("Unknown variant caller for merging SBSs into DBSs",
         variant.caller,
         "\nDid you want to use Read_DBS_VCF_and_BAMs_to_verify_DBSs")
  }

  if (mode(input.vcf) == "character" && is.null(dim(input.vcf))) {
    # Interpret this as a file path
    input.vcf.name <- input.vcf

    if (!file.exists(input.vcf.name)) {
      stop("VCF file ", input.vcf.name, "does not exist")
    }

    vcf.list <-ICAMS::ReadAndSplitVCFs(input.vcf.name,
                                       variant.caller   = variant.caller,
                                       num.of.cores     = 1,
                                       get.vaf.function = NULL, # get.vaf.function,
                                       max.vaf.diff     = 1,
                                       always.merge.SBS = TRUE)

    input.vcf <- vcf.list$DBS[[1]]

    if (verbose > 0) {
      message("Analyzing vcf            ", input.vcf.name)
    }

    if (is.null(outfile)) {
      outfile <- paste0(input.vcf.name, "_evaluated.vcf")
    }
  } else {
    stop("Argument input.vcf to Read_SBS_VCF_and_BAMs_to_verify_DBSs ",
         "must be a file path")
  }

  eval.out <-
    Read_DBS_VCF_and_BAMs_to_verify_DBSs(input.vcf        = input.vcf,
                                         Nbam.name        = Nbam.name,
                                         Tbam.name        = Tbam.name,
                                         N.slice.dir      = N.slice.dir,
                                         T.slice.dir      = T.slice.dir,
                                         unlink.slice.dir = unlink.slice.dir,
                                         verbose          = verbose,
                                         outfile          = outfile)
  invisible(eval.out)
}
