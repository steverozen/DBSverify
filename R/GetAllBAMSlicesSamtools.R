#' Get and optionally save all slices of BAM as specified by positions in a VCF table.
#'
#' @param vcf A VCF (Variant Call Format) "file" as a data.table or similar.
#'
#' @param bam.name The file name (path) to the BAM file to slice.
#'
#' @param padding How many base pairs on each side of the first base of the DBS
#'    to to keep in the BAM slices.
#'
#' @param where.to.put.slices If \code{NULL}, create a temporary directory to store
#'    the slices.  Otherwise, a character string that specifies a directory in
#'    which to store the BAM slices. This is directory is created if necessary.
#'
#' @return A character string that specifies the directory the containing the
#'   BAM slices. The slices are stored as SAM files.

GetAllBAMSlicesSamtools <- function(vcf,
                                    bam.name,
                                    padding             = 10,
                                    where.to.put.slices = tempfile()) {
  if (!dir.exists(where.to.put.slices)) {
    if (file.exists(where.to.put.slices)) {
      stop(where.to.put.slices, " exists as a file; needs to be a directory")
    }
    if (!dir.create(where.to.put.slices, recursive = TRUE)) {
      stop("Unable to create directory ", where.to.put.slices)
    }
  }
  cat("This directory contains BAM slices from",
      bam.name, "\n", as.character(Sys.time()),"\n",
      file = file.path(where.to.put.slices, "README.txt"))

  for(v in 1:nrow(vcf)){

    CHROM <- vcf[v, "CHROM"]
    POS   <- vcf[v, "POS"]
    REF   <- vcf[v, "REF"]
    ALT   <- vcf[v, "ALT"]

    slice.filename <- file.path(where.to.put.slices,
                                paste0(CHROM, "-", POS, ".sam"))
    CreateBAMSliceFileSamtools(bam.name,
                               CHROM,
                               POS,
                               padding = padding,
                               save.file.path = slice.filename)
  }
  return(where.to.put.slices)
}
