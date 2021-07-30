#' Check whether analysis of individual reads supports DBSs in a VCF.
#'
#' @param vcf An in-memory representation of a "variant call file", VCF, as a \code{data.frame}.
#'
#' @param Nbam.name The name of the BAM file for the normal sample.
#'
#' @param Tbam.name The name of the BAM file for the tumor sample.
#'
#' @return An in-memory VCF based on the input \code{vcf} with the additional columns
#'   \code{ID}, \code{NreadSupport}, \code{TreadSupport}, and \code{DBSconclusion}.
#'   If \code{vcf} is empty, return it without adding columns.
#'
#'   XXXXXX MORE MORE
#'
#' @export

VerifyDBSVcf_old <- function(vcf, Nbam.name, Tbam.name) {

  CheckBAM(Nbam.name)
  CheckBAM(Tbam.name)

  if (nrow(vcf) == 0) return(vcf)

  # debug(addID)
  vcf$ID <- IDVector(vcf)
  ## remove duplicates if present
  vcf<-vcf[!duplicated(vcf$ID), ]
  rownames(vcf)<-vcf$ID

  vcf$Format<-paste0("WtReads:pos1reads:pos2reads:MutReads")

  vcf2 <- SummarizeReadSupport(
    vcf = vcf, Nbam.name = Nbam.name, Tbam.name = Tbam.name)
  #use system2

  vcf3 <-DBSConclusion(vcf2)
  return(vcf3)
}
