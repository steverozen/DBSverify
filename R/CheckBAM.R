CheckBAM <- function(bam.name) {
  if (!file.exists(bam.name)) {
    stop("BAM file ", bam.name, " does not exist" )
  }

  index <- paste0(bam.name, ".bai")
  if (!file.exists(index)) {
    stop("BAM index file ", index, " does not exist; run samtools index" )
  }

}
