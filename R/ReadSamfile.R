#' Read a "SAM file", discarding some reads that cannot be interpreted for our purposes
#'
#' @details
#' SAM stands for "Sequence Alignment Map", a text file that
#' represents aligned next-generation sequencing reads
#' (https://en.wikipedia.org/wiki/SAM_(file_format)).
#' Only keep reads that statisfy certain conditions:
#'
#' * Mate pair maps to same chromosome
#' * Mapping quality >= 30
#' * "FLAG" < 256
#'
#' @param filename The name of the SAM file to read.
#'
#' @param check.CIGAR Include only reads for which the CIGAR string
#'  is only \code{\\d+M} (one or more digits followed M, and nothing else before or after).
#'  This means there are no insertions or deleitons in the read
#'  versus the reference and there is no soft clipping.
#'
#' @return A \code{data.frame} with \code{colnames} for the first 11
#' columns, with \code{colnames} as specified in
#' https://en.wikipedia.org/wiki/SAM_(file_format)
#' and one row per read.
#'
#' @keywords internal

ReadSamfile <- function(filename, check.CIGAR = TRUE) {
  df <- readLines(filename)

  # SAM file headers start with "@"; discard them
  df <- df[-grep("^@", df)]

  df <- as.data.frame(stringr::str_split_fixed(df, "\t", 16))
  colnames(df)[1:11] <-
    c("QNAME", "FLAG", "CHROM", "POS", "MAPQ", "CIGAR", "Mate_CHROM",
      "Mate_POS", "InsertSize", "SEQ", "QUAL")

  df.n <- nrow(df)
  if (Sys.getenv("TRACEDBS") != "") {
    message("\nReadSamfile: filename = ", filename)
    message("ReadSamfile: starting with ", df.n, " reads")
  }

  # remove weird rows if required
  df <- df[!is.na(suppressWarnings(as.numeric(df$FLAG))), ]
  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: ",
            df.n - nrow(df),
            " reads removed because FLAG could not interpreted as a number")
  }
  df.n <- nrow(df)

  # Include only reads with FLAG < 256
  # Higher values mark reads that
  # - failed vendor QC
  # - have multiple alignments
  # - are marked as duplicates
  bad.flag <- which(as.numeric(df$FLAG) >= 256)
  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: ",
            length(bad.flag),
            " reads with bad FLAGs removed")
    if (length(bad.flag) > 0) {
      message("bad FLAGs: ",
              paste(df[bad.flag, "FLAG"], collapse = ", "))
    }
  }
  df <- df[as.numeric(df$FLAG) < 256, ]

  # Keep only reads that have a pair on the same chromosome.
  df <- df[df$Mate_CHROM == "=", ]
  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: ",
            df.n - nrow(df),
            " reads removed pair was on another chromosome")
  }
  df.n <- nrow(df)

  df <- df[df$MAPQ >= 30, ]
  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: ",
            df.n - nrow(df),
            " reads removed because MAPQ < 30")
  }
  df.n <- nrow(df)

  if (check.CIGAR) {
    cigar.good <- grep("^\\d+M$", df$CIGAR)
    cigar.bad <- df[-cigar.good, ]
    df <- df[cigar.good, ]
    if (Sys.getenv("TRACEDBS") != "") {
      message("ReadSamfile: ",
              df.n - nrow(df),
              " reads removed because CIGAR string was not \\d+M")
      message("Bad CIGARs: ", paste(cigar.bad$CIGAR, collapse = ", "))
    }
  }
  df.n <- nrow(df)
  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: returning ", df.n, " reads")
  }

  return(df)
}
