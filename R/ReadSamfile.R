#' Read a "SAM file", discarding some reads that cannot be interpreted for our purposes
#'
#' SAM stands for "Sequence Alignment Map", a text file that
#' represents aligned next-generation sequencing reads
#' (https://en.wikipedia.org/wiki/SAM_(file_format))
#'
#' @param filename The name of the SAM file to read.
#'
#' @return A \code{data.frame} with \code{colnames} for the first 11
#' columns, with \code{colnames} as specified in
#' https://en.wikipedia.org/wiki/SAM_(file_format)
#' and one row per read.
#'
#' @keywords internal

ReadSamfile <- function(filename) {
  df <- utils::read.csv(
    filename, sep = "^", as.is = TRUE, fill = TRUE, header = FALSE)

  # SAM file headers start with "@"; discard them
  df <- df[-grep("^@", df[, 1]), ]

  df <- as.data.frame(stringr::str_split_fixed(df, "\t", 16))
  colnames(df)[1:11] <-
    c("QNAME", "FLAG", "CHROM", "POS", "MAPQ", "CIGAR", "Mate_CHROM",
      "Mate_POS", "InsertSize", "SEQ", "QUAL")

  ## remove weird rows if required
  df <- df[!is.na(suppressWarnings(as.numeric(df$FLAG))), ]

  ## remove reads with flag > 256
  ## these are reads that
  ## - fail vendor QC
  ## - have multiple alignments
  ## - are marked as duplicates
  ## - or their pair aligns to a different chromosome
  df <- df[as.numeric(df$FLAG) < 256, ]
  df <- df[df$Mate_CHROM == "=", ]
  return(df)
}
