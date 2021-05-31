#' Get a BAM slice as a \code{data.frame} possibly using different functions.
#'
#' @param BAM.name The input BAM file name.
#'
#' @param CHROM The chromosome from which to slice.
#'
#' @param POS The position in \code{CHROM} from which to slice.
#'
#' @param padding The number of basepairs on each side of \code{POS} to take.
#'
#' @param read.fn The function to use to get the BAM slice as a SAM file.
#'
#' @param ... Additional arguments to read.fn
#'
#' @return The SAM file contents as a \code{data.frame}.

CreateBAMSlice <- function(BAM.name, CHROM, POS, padding = 10, read.fn = CreateBAMSliceFileSamtools, ...) {
  sam.file <- read.fn(BAM.name, CHROM, POS, padding, ...)
  sam <- ReadSamfile(sam.file)
  if (0 != unlink(sam.file)) warning("unable to unlink ", sam.file)
  return(sam)
}
