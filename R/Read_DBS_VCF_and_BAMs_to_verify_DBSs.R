#' Determine whether sequencing reads in fact support DBSs present in a VCF file.
#'
#' @param input.vcf If a character string, then the path to a VCF file; otherwise
#'   A a single VCF "file" as a data.frame or similar object.
#'
#' @param Nbam.name The name of the BAM file for the normal sample corresponding to \code{vcf.name}.
#'
#' @param Tbam.name The name of the BAM file for the tumor sample corresponding to \code{vcf.name}.
#'
#' @param N.slice.dir Directory for the slices of the normal BAM.
#'  Created if necessary.
#'
#' @param T.slice.dir Directory for the slices of the tumor BAM.
#'  Created if necessary. Must be different than \code{N.slice.dir}.
#'
#' @param unlink.slice.dir If \code{TRUE} unlink \code{N.slice.dir}
#'  and \code{T.slice.dir} before return.
#'
#' @param exclude.SBSs If \code{TRUE} silently filter out (exclude)
#'    SBSs in the input VCF. This makes sense if the the VCF is from
#'    a caller (like Mutect or the Hartwig Medical Foundation caller)
#'    that calls both SBSs and DBS.
#'
#' @param verbose If > 0 print a message when starting the number of slices
#'    generated every \code{verbose} slices.
#'
#' @param outfile If not \code{NULL} then write the "evaluated" VCF to \code{outfile};
#'  otherwise write it to \code{paste0(input.vcf(vcf.name, "_evaluated.vcf")}. Must be
#'  non-\code{NULL} if \code{} is not a file path.
#'
#' @details Creates a new VCF file.
#'  This VCF file has no data rows if there were no DBSs to analyze.
#'  Otherwise, this VCF contains the additional columns.
#'  Any SBSs or indels in the input are silently ignored,
#'  and no attempt is made to merge adjacent SBSs.
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
#'     believable (\code{"True DBS"}), or if the DBS is not
#'     believable, a string that describes
#'     why not.
#'
#' @return Invisibly, a list with the elements
#'
#' 1. The name of the DBS-only VCF file created.
#'
#' 1. The in-memory representation of the DBS VCF as a \code{data.table}.
#'
#' 1. The name of the directory with the normal sam slices, if \code{unlink.slice.dir} is \code{FALSE}.
#'
#' 1. The name of the directory with the tumor sam slices, if \code{unlink.slice.dir} is \code{FALSE}.
#'
#' @md
#' @export
#'

Read_DBS_VCF_and_BAMs_to_verify_DBSs <- function(input.vcf,
                                                 Nbam.name,
                                                 Tbam.name,
                                                 N.slice.dir      = tempfile(),
                                                 T.slice.dir      = tempfile(),
                                                 unlink.slice.dir = TRUE,
                                                 exclude.SBSs     = TRUE,
                                                 verbose          = 0,
                                                 outfile          = NULL) {

  CheckBAM(Nbam.name)
  CheckBAM(Tbam.name)
  # TestBAMAndSamtools(Nbam.name) # generates spurious warnings when working with minibams
  # TestBAMAndSamtools(Tbam.name)
  if (N.slice.dir == T.slice.dir) {
    stop("T.slice.dir and N.slice.dir must be different; got ",
         N.slice.dir)
  }

  if (mode(input.vcf) == "character" && is.null(dim(input.vcf))) {
    # Interpret this as a file path
    input.vcf.name <- input.vcf
    vcf.list <-ICAMS::ReadAndSplitVCFs(input.vcf.name,
                                       variant.caller   = "unknown",
                                       num.of.cores     = 1,
                                       alwaysa.merge.SBS = FALSE)

    input.vcf <- vcf.list$DBS[[1]]

    if (verbose > 0) {
      message("Analyzing vcf            ", input.vcf.name)
    }

    if (is.null(outfile)) {
      outfile <- paste0(input.vcf.name, "_evaluated.vcf")
    }
  } else {
    if (is.null(ncol(input.vcf))) {
      stop("The argument input.vcf does not have any columns; looks like a scalar")
    }
    if (is.array(input.vcf)) {
      input.vcf <- data.frame(input.vcf)
    }
    if (is.null(outfile)) {
      stop("Please provide a value for argument outfile")
    }
  }

if (verbose > 0) {
    message("Normal BAM               ", Nbam.name)
    message("Tumor BAM                ", Tbam.name)
  }

  input.vcf <- Remove_non_canonical_chromosomes(input.vcf, verbose)
  if (nrow(input.vcf) == 0) { message("No DBSs found") }

  evaluated.vcf <- VerifyDBSVcf(input.vcf,
                                Nbam.name = Nbam.name,
                                Tbam.name = Tbam.name,
                                N.slice.dir = N.slice.dir,
                                T.slice.dir = T.slice.dir,
                                verbose     = verbose)

 if (unlink.slice.dir) {
    unlink(N.slice.dir, recursive = TRUE)
    unlink(T.slice.dir, recursive = TRUE)
    N.slice.dir <- NA
    T.slice.dir <- NA
  }

  message("Output evaluated VCF file ", outfile)

  cat("#", file = outfile)
  suppressWarnings(
    utils::write.table(
      evaluated.vcf,
      file = outfile,
      sep="\t",
      quote = FALSE,
      row.names = FALSE,
      append = TRUE))

  invisible(list(evaluated.vcf.name  = outfile,
                 evaluated.vcf       = evaluated.vcf,
                 N.slice.dir         = N.slice.dir,
                 T.slice.dir         = T.slice.dir))
}
