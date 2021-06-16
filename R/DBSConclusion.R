#' Examines every row of DBS VCF already processed by \code{\link{Slice2ReadSupport}} and decides which DBSs are real.
#'
#' @param vcf An in memory DBS VCF (as a \code{data.frame})
#' already processed by \code{\link{Slice2ReadSupport}}  (so that the fields
#' \code{vcf$NreadSupport} and \code{vcf$TreadSupport} are populated).
#'
#' @param germlineCutOff If this proportion of normal reads show one or the
#'   other variant (or both variants), consider this a germline variant or
#'   partial germline variant.
#'
#' @param max.non.mut.reads Tolerate this number of reads in the tumor
#'   that do not support both mutated positions.
#'
#' @return A vcf with the field \code{DBSconclusion} populated.
#'
#' @keywords internal

DBSConclusion <- function(vcf, germlineCutOff = 0.2, max.non.mut.reads = 1) {
  if (nrow(vcf) == 0) return(vcf)

  rr <- apply(X = vcf,
              MARGIN = 1,
              FUN = DBS_conclusion_1_row,
              germlineCutOff = germlineCutOff,
              max.non.mut.reads = max.non.mut.reads)

  vcf$DBSconclusion <- unlist(rr)

  return(vcf)
}
