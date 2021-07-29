#' Examine 1 row of DBS VCF already processed by \code{\link{Slice2ReadSupport}} and decide which DBSs are real.
#'
#' @param row One row of a DBS VCF (as a \code{data.frame})
#' already processed by \code{\link{Slice2ReadSupport}}  (so that the fields
#' \code{vcf$NreadSupport} and \code{vcf$TreadSupport} are populated).
#'
#' @param germlineCutOff If this proportion of normal reads show one or the
#'   other variant (or both variants), consider this a germline variant or
#'   partial germline variant.
#'
#' @param max.half.support.T.reads Do not tolerate more than this number of reads in the tumor
#'   that support one but not both mutated positions.
#'
#' @return A character string indicating the conclusion about the putative DBS.

DBS_conclusion_1_row <- function(row, germlineCutOff = 0.2, max.half.support.T.reads = 1) {
  wt   <- 1
  pos1 <- 2
  pos2 <- 3
  dbs  <- 4
  p.cutoff <- 0.05
  N.read.counts <-
    as.numeric(unlist(strsplit(row["NreadSupport"], ":", fixed = TRUE)))

  xx <- sum(N.read.counts)
  if (xx == 0) {
    return("No high quality normal reads")
  }

  N.read.prop   <- N.read.counts / xx

  T.read.counts <-
    as.numeric(unlist(strsplit(row["TreadSupport"], ":", fixed = TRUE)))
  #if (N.read.prop[pos1] >= germlineCutOff &&
  #    N.read.prop[pos2] >= germlineCutOff) {
  #  return("DBS overlaps germline SNPs")
  # }
  if (N.read.prop[pos1] >= germlineCutOff) {
    return("DBS overlaps germline SNP")
  }
  if (N.read.prop[pos2] >= germlineCutOff) {
    return("DBS overlaps germline SNP")
  }

  if (N.read.prop[dbs] >= germlineCutOff) {
    return("Germline DBS")
  }

  if (N.read.counts[dbs] > 1) {
    return("> 1 normal read with DBS")
  }

  if (xx < 5) {
    return("< 5 high quality normal reads")
  }

  if (sum(T.read.counts[c(pos1, pos2, dbs)]) == 0) {
    return("Neither position supported")
  }

  T.dbs <- T.read.counts[dbs]
  if (sum(T.read.counts[c(pos1, pos2)]) > max.half.support.T.reads) {
    return("Adjacent SBSs")
  } else {
    if (T.dbs == 0) {
      return("0 tumor reads support the DBS")
    }
    if (T.dbs == 1) {
      return("Only 1 tumor read supports the DBS")
    }
    # At this point T.dbs > 1) and N.read.counts[dbs] %in% 0:1
    if (as.numeric(row["num_bad_mapped_reads"]) >= sum(T.read.counts)) {
      return("Too many badly mapped tumor reads")
    }
    if (as.numeric(row["num_bad_mapped_DBS_reads"]) >= T.dbs) {
      return("Too many badly mapped DBS reads")
    }

    if (N.read.counts[dbs] == 0) {
      return("True DBS")
    }
    if (N.read.counts[dbs] == 1) {
      # Could the DBS reads just be the same noise
      # as in the normal?
      pp <- stats::fisher.test(
        matrix(c(N.read.counts[wt], 1,
                 T.read.counts[wt], T.dbs),
               ncol = 2),
        alt = "g")$p.value
      if (pp > p.cutoff) {
        return("Failed Fisher test")
      } else {
        return("True DBS")
      }
    }
  }

  return(paste("Should not get here",
               row["NreadSupport"],
               row["TreadSupport"]))
}
