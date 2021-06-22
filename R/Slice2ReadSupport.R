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

  sam.df <- ReadSamfile(file.path(slice.dir,paste0(CHROM, "-", POS, ".sam")))

  if (nrow(sam.df) == 0) {
    return("0:0:0:0")
  }

  categories <- CategorizeReads(sam = sam.df, POS = POS, REF = REF, ALT = ALT)

  read.support <-paste0(sum(categories == "WT read"),":",
                        sum(categories == "Read supports only 1st position"),":",
                        sum(categories == "Read supports only 2nd position"),":",
                        sum(categories == "Mut read"))

  return(read.support)

}
