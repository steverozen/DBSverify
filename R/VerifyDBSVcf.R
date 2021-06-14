#' Check whether analysis of individual reads supports DBSs in a VCF.
#'
#' @param vcf An in-memory representation of a "variant call file", VCF, as a \code{data.frame}.
#'
#' @param Nbam.name The name of the BAM file for the normal sample.
#'
#' @param Tbam.name The name of the BAM file for the tumor sample.
#'
#' @param N.slice.dir The directory containing the slices of the normal BAM.
#'
#' @param T.slice.dir The directory containing the slices of the tumor BAM. Must
#'    be different than \code{N.slice.dir}.
#'
#' @param verbose If > 0 print a message when starting the number of slices
#'    generated every \code{verbose} slices.
#'
#' @param padding The number of base pairs on either side of the first position
#'  of the DBS to include
#'  in the slices.
#'
#' @return An in-memory VCF based on the input \code{vcf} with the additional columns
#'   as described in \code{\link{Read_DBS_VCF_and_BAMs_to_verify_DBSs}}.
#'
#' @keywords internal

VerifyDBSVcf <- function(vcf,
                         Nbam.name,
                         Tbam.name,
                         N.slice.dir,
                         T.slice.dir,
                         padding = 10,
                         verbose = 0) {

  CheckBAM(Nbam.name)
  CheckBAM(Tbam.name)
  if (N.slice.dir == T.slice.dir) {
    stop("T.slice.dir and N.slice.dir must be different; got ",
         N.slice.dir)
  }

  if (nrow(vcf) == 0) return(vcf)

  vcf$ID <- IDVector(vcf)
  vcf<-vcf[!duplicated(vcf$ID), ]
  rownames(vcf)<-vcf$ID

  vcf$Format<-paste0("WtReads:pos1reads:pos2reads:MutReads")

  GetAllBAMSlicesSamtools(vcf                 = vcf,
                          bam.name            = Nbam.name,
                          padding             = padding,
                          where.to.put.slices = N.slice.dir,
                          verbose             = verbose)

  GetAllBAMSlicesSamtools(vcf                 = vcf,
                          bam.name            = Tbam.name,
                          padding             = padding,
                          where.to.put.slices = T.slice.dir,
                          verbose             = verbose)

  vcf2 <- SummarizeReadSupportFromSlices(vcf, N.slice.dir, T.slice.dir)

  vcf3 <-DBSConclusion(vcf2)

  return(vcf3)
}
