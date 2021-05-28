#' Create a BAM slice as a SAM file using the samtools executable.
#'
#' @param BAM.name The input BAM file name.
#'
#' @param CHROM The chromosome from which to slice.
#'
#' @param POS The position in \code{CHROM} from which to slice.
#'
#' @param padding The number of basepairs on each side of \code{POS} to take.
#'
#' @param save.file.path If not NULL, save the file to that file path.
#'
#' @return The path of the SAM file that contains the BAM slice.

CreateBAMSliceFileSamtools <-
  function(BAM.name, CHROM, POS, padding = 10, save.file.path = NULL) {
    BAM.coord <- paste0(CHROM, ":", POS - padding, "-", POS + padding)
    if (is.null(save.file.path)) {
      save.file.path <- tempfile(pattern = "")
    }
    system(paste0("samtools view -h ", BAM.name, " ",
                  BAM.coord, " > ", save.file.path), wait = T)
    return(save.file.path)
  }
