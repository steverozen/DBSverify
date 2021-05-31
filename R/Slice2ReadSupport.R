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

Slice2ReadSupport <- function(slice.dir, CHROM, POS, REF, ALT) {

  sam <- ReadSamfile(file.path(slice.dir,paste0(CHROM, "-", POS, ".sam")))

  eval <- CategorizeReads(sam = sam, POS = POS, REF = REF, ALT = ALT)

  read.support <-paste0(sum(eval == "WT read"),":",
                        sum(eval == "Read supports only 1st position"),":",
                        sum(eval == "Read supports only 2nd position"),":",
                        sum(eval == "Mut read"))

  return(read.support)

}
