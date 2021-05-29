#' Calculate the read support for putative DBS a given position in a BAM file
#'
#' @inheritParams CreateBAMSlice
#'
#' @param REF The reference variant
#'
#' @param ALT The alternate variant

GetReadSupport <- function(BAM.name, CHROM, POS, REF, ALT,
                           read.fn = CreateBAMSliceFileSamtools) {

  sam <- CreateBAMSlice(BAM.name = BAM.name, CHROM = CHROM, POS = POS, read.fn = read.fn)

  eval <- CategorizeReads(sam = sam, POS = POS, REF = REF, ALT = ALT)

  read.support <-paste0(sum(eval == "WT read"),":",
                        sum(eval == "Read supports only 1st position"),":",
                        sum(eval == "Read supports only 2nd position"),":",
                        sum(eval == "Mut read"))

  return(read.support)

}
