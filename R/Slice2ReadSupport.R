#' Calculate the support for a putative DBS from the SAM slice containing the overlapping reads.
#'
#' @param slice.dir The directory containing the sam / bam slices.
#'
#' @param CHROM The chromosome identifier.
#'
#' @param POS The first position of the DBS.
#'
#' @param REF The reference variant.
#'
#' @param ALT The alternate variant.
#'
#' @keywords internal

Slice2ReadSupport <- function(slice.dir, CHROM, POS, REF, ALT) {

  stopifnot(!is.na(REF))
  stopifnot(!is.na(ALT))

  sorted.reads <- ReadSamfile(file.path(slice.dir,paste0(CHROM, "-", POS, ".sam")))

  good.reads <- sorted.reads$good.reads

  bad.mapped.reads <- rbind(sorted.reads$reads.with.bad.MAPQ,
                            sorted.reads$reads.with.bad.Mate_CHROM)

  if (nrow(bad.mapped.reads) > 0 ) {
    bad.categories <-
      CategorizeReads(sam = bad.mapped.reads, POS = POS, REF = REF, ALT = ALT)
    num.bad.mapped.DBS.reads <- sum(bad.categories == "Mut read")
  } else {
    num.bad.mapped.DBS.reads <- 0
  }

  if (nrow(good.reads) == 0) {
    return(list(read.support             = "0:0:0:0",
                num.bad.mapped.reads     = nrow(bad.mapped.reads),
                num.bad.mapped.DBS.reads = num.bad.mapped.DBS.reads))
  }

  categories <- CategorizeReads(sam = good.reads, POS = POS, REF = REF, ALT = ALT)

  read.support <-paste0(sum(categories == "WT read"),":",
                        sum(categories == "Read supports only 1st position"),":",
                        sum(categories == "Read supports only 2nd position"),":",
                        sum(categories == "Mut read"))

  return(list(read.support             = read.support,
              num.bad.mapped.reads     = nrow(bad.mapped.reads),
              num.bad.mapped.DBS.reads = num.bad.mapped.DBS.reads))

}
