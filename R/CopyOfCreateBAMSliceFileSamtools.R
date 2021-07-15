#' Create a BAM slice as a SAM file using the samtools executable.
#'
#' @param BAM.name The input BAM file name.
#'
#' @param CHROM The chromosome from which to slice.
#'
#' @param POS The position in \code{CHROM} from which to slice.
#'
#' @param padding The number of base pairs on each side of \code{POS} to take.
#'
#' @param save.file.path Where to save the sam file.
#'
#' @return The path of the SAM file that contains the BAM slice.
#'
#' @keywords internal

CreateBAMSliceFileSamtools <-
  function(BAM.name, CHROM, POS, padding = 10, save.file.path) {
    POS <- as.numeric(POS) # In case it was a character string.
    BAM.coord <- paste0(CHROM, ":", POS - padding, "-", POS + padding)
    if (is.null(save.file.path)) {
      save.file.path <- tempfile(pattern = "")
    }
    status <- system2("samtools",
                      c("view", "-h", BAM.name, BAM.coord),
                      stdout = TRUE,
                      wait = TRUE,
                      stderr = TRUE)
    status2 <- attr(status, "status")
    errmsg <- attr(status, "errmsg")
    if (!is.null(errmsg)) {
      stop("samtools generated errmsg ", errmsg)
    }
    if (!is.null(status2)) {
      if (status2 == 127) stop("samtools: command not found")
      if (status2 != 0) stop("samtools returned error status ", status2)
    }

    if (!grepl(pattern = "^@", status[1])) {
      message("POSSIBE PROBLEM, first line of sam slice does not start with @")
    }
    cat(status, sep = "\n", file = save.file.path)

    return(save.file.path)
  }
