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

old_DBSConclusion <- function(vcf, germlineCutOff = 0.2, max.non.mut.reads = 1) {
  if (nrow(vcf) == 0) return(vcf)

  vcf$DBSconclusion <- NA

  x.data.frame <- function(read.support) {
    tmp <- as.data.frame(stringr::str_split_fixed(read.support, ":", 4))
    for(cc in 1:4) tmp[ , cc] <- as.numeric(tmp[ , cc])
    # tmp <- apply(tmp, MARGIN = 2, FUN = as.numeric) # apply
    colnames(tmp) <- c("WT", "pos1only", "pos2only", "mutant")
    # tmp isa data frame of counts of reads of each category, WT, pos1only, ...
    return(tmp)
  }

  norm <- x.data.frame(vcf$NreadSupport)

  i <- sweep(norm, 1, rowSums(norm), FUN = "/") # i is proportions of total reads

  vcf$DBSconclusion[i[, "pos1only"] >= germlineCutOff & i[, "pos2only"] >= germlineCutOff] <-
    "DBS overlaps 2 germline SNPs"

  vcf$DBSconclusion[i[, "pos1only"] >= germlineCutOff & i[, "pos2only"] < germlineCutOff] <-
    "1st position overlaps germline SNP"

  vcf$DBSconclusion[i[, "pos1only"] < germlineCutOff & i[, "pos2only"] >= germlineCutOff] <-
    "2nd position overlaps germline SNP"

  vcf$DBSconclusion[i[, "mutant"] >= germlineCutOff] <- "Germline DBS"

  # Whatever is NA does not overlap a germline SNP
  i <- is.na(vcf$DBSconclusion)

  mut <- x.data.frame(vcf$TreadSupport)

  vcf$DBSconclusion[i & rowSums(mut[, 2:3]) <= max.non.mut.reads &
                      mut[, "mutant"] > 0] <-
    "True DBS"

  vcf$DBSconclusion[i & rowSums(mut[, 2:4]) == 0] <- "Neither position supported"

  vcf$DBSconclusion[i & rowSums(mut[, 2:3]) > 1] <- "Adjacent SBSs"

  return(vcf)
}
