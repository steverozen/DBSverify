#' @keywords internal

CheckBAM <- function(bam.name, must.succeed = TRUE) {
  if (!file.exists(bam.name)) {
    if (must.succeed) {
      stop("BAM file ", bam.name, " does not exist" )
    } else {
      return(FALSE)
    }
  }

  index <- paste0(bam.name, ".bai")
  if (!file.exists(index)) {
    stop("BAM index file ", index, " does not exist; run samtools index" )
  }

  return(TRUE)
}
