#' Agressively merge SBSs into DBSs then determine whether sequencing reads in fact support DBSs.
#'
#' @param variant.caller One of \code{"strelka"}, \code{"PCAWG"},
#'  or \code{"unknown"}.
#'  Merging adjacent SBS is done by \code{\link[ICAMS]{ReadAndSplitVCFs}}.
#'  Do not use this for callers that
#'
#' @inheritParams Read_DBS_VCF_and_BAMs_to_verfiy_DBSs
#'
#' @details Creates a new VCF file.
#'  This VCF file has no data rows if there were no DBSs to analyze.
#'  Otherwise, this VCF contains the additional columns
#'
#'  1. \code{Format} contains the fixed string \code{"WtReads:pos1reads:pos2reads:MutReads"}.
#'
#'  1. \code{NreadSupport} With regard to the two positions of the DBS in
#'     the normal BAM, a string with 4 numbers separated by ":", with the numbers
#'     indicating respectively:
#'
#'     * the number of reads that are reference sequence at
#'     both positions of the DBS,
#'
#'     * the number of reads that that have the alternative
#'     allele only at the 1st position of the DBS,
#'
#'     * the number of reads that
#'     have the alternative allele only at the second position of the DBS, and
#'
#'     * the number
#'     of reads that have the alternative alleles at both positions of the DBS.
#'
#'  1. \code{TreadSupport} Information analogous to that in \code{NreadSupport}, for the
#'     tumor BAM.
#'
#'  1. \code{DBSconclusion} A string that describes whether the DBSs is
#'     believable (\code{"True DBS"}), and if not, a string that describes
#'     why not.
#'
#' @return Invisibly, a list with the elements
#'
#' 1. The name of the DBS-only VCF file created.
#'
#' 1. The in-memory representation of the DBS VCF as a \code{data.table}.
#'
#' 1. The directory with the normal sam slices, if \code{unlink.slice.dir} is \code{FALSE}.
#'
#' 1. The directory with the tumor sam slices, if \code{unlink.slice.dir} is \code{FALSE}.
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
  if (N.slice.dir == T.slice.dir) {
    stop("T.slice.dir and N.slice.dir must be different; got ",
         N.slice.dir)
  }

  if (variant.caller %in% c("strelka", "PCAWG", "unknown")) {
    get.vaf.function <- NULL
    if (variant.caller == "PCAWG") {
      variant.caller <- "unknown"
      get.vaf.function <- GetPCAWGVAF
    }
  } else {
    stop("Unknown variant caller for merging SBSs into DBSs",
         variant.caller,
         "\nDid you want to use Read_DBS_VCF_and_BAMs_to_verfiy_DBSs")
  }

  if (mode(input.vcf) == "character") {
    # Interpret this as a file path
    input.vcf.name <- input.vcf

    if (!file.exists(input.vcf.name)) {
      stop("VCF file ", input.vcf.name, "does not exist")
    }

    vcf.list <-ICAMS::ReadAndSplitVCFs(input.vcf.name,
                                       variant.caller   = variant.caller,
                                       num.of.cores     = 1,
                                       get.vaf.function = get.vaf.function,
                                       max.vaf.diff     = 1,
                                       always.merge.SBS = TRUE)

    input.vcf <- vcf.list$DBS[[1]]

    if (verbose > 0) {
      message("Analyzing vcf            ", input.vcf.name)
    }

    if (is.null(outfile)) {
      outfile <- paste0(input.vcf.name, "_evaluated.vcf")
    }
  }

  eval.out <-
    Read_DBS_VCF_and_BAMs_to_verfiy_DBSs(input.vcf        = input.vcf,
                                         Nbam.name        = Nbam.name,
                                         Tbam.name        = Tbam.name,
                                         N.slice.dir      = N.slice.dir,
                                         T.slice.dir      = T.slice.dir,
                                         unlink.slice.dir = unlink.slice.dir,
                                         verbose          = verbose,
                                         outfile          = outfile)
  invisible(eval.out)
}
