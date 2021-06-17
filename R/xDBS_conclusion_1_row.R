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
#' @param max.non.mut.reads Tolerate this number of reads in the tumor
#'   that do not support both mutated positions.
#'
#' @return A character string indicating the conclusion about the putative DBS.

xDBS_conclusion_1_row <- function(row, germlineCutOff = 0.2, max.non.mut.reads = 1) {
  wt   <- 1
  pos1 <- 2
  pos2 <- 3
  dbs  <- 4
  p.cutoff <- 0.05
  N.read.counts <-
    as.numeric(unlist(strsplit(row["NreadSupport"], ":", fixed = TRUE)))

  ss <- sum(N.read.counts)
  if (ss == 0) {
    return("No high quality normal reads")
  }

  N.read.prop   <- N.read.counts / ss

  T.read.counts <-
    as.numeric(unlist(strsplit(row["TreadSupport"], ":", fixed = TRUE)))
  if (N.read.prop[pos1] >= germlineCutOff &&
      N.read.prop[pos2] >= germlineCutOff) {
    return("DBS overlaps 2 germline SNPs")
  }
  if (N.read.prop[pos1] >= germlineCutOff) {
    return("1st position overlaps germline SNP")
  }
  if (N.read.prop[pos2] >= germlineCutOff) {
    return("2nd position overlaps germline SNP")
  }

  if (N.read.prop[dbs] >= germlineCutOff) {
    return("Germline DBS")
  }

  if (N.read.counts[dbs] > 1) {
    return("> 1 normal read with DBS")
  }

  if (ss < 5) {
    return("< 5 high quality normal reads")
  }

  if (sum(T.read.counts[c(pos1, pos2, dbs)]) == 0) {
    return("Neither position supported")
  }

  # It could be e.g. 0:1:0 for pos1 pos2 dbs

  if ((sum(T.read.counts[c(pos1,pos2)]) <= max.non.mut.reads)) {
    if (T.read.counts[dbs] == 0) {
      return("0 tumor reads support the DBS")
    }
    if (T.read.counts[dbs] == 1) {
      return("Only 1 tumor read supports the DBS")
    }
    if (N.read.counts[dbs] == 0) {
      return("True DBS")
    }
    if (N.read.counts[dbs] == 1) {
      pp <- stats::fisher.test(
        matrix(c(N.read.counts[wt], 1,
                 T.read.counts[wt], T.read.counts[dbs]),
               ncol = 2),
        alt = "g")$p.value
      if (pp > p.cutoff) {
        return("Proporiton of tumor reads with DBS too low")
      } else {
        return("True DBS")
      }
    }
  }

   if (sum(T.read.counts[c(pos1, pos2)]) > 1) { # 1 should be max.non.mut.reads
    return("Adjacent SBSs")
  }

  return(paste("Should not get here",
               row["NreadSupport"],
               row["TreadSupport"]))


}
