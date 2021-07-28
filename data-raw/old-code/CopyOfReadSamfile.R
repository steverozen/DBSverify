#' Read a "SAM file", discarding some reads that cannot be interpreted for our purposes
#'
#' @details
#' SAM stands for "Sequence Alignment Map", a text file that
#' represents aligned next-generation sequencing reads
#' (https://en.wikipedia.org/wiki/SAM_(file_format)).
#' Only keep reads that satisfy certain conditions:
#'
#' * Mate pair maps to same chromosome
#' * Mapping quality >= 30
#' * "FLAG" < 256
#'
#' @param filename The name of the SAM file to read.
#'
#' @param check.CIGAR Include only reads for which the CIGAR string
#'  is only \code{\\d+M} (one or more digits followed M, and nothing else before or after).
#'  This means there are no insertions or deletions in the read
#'  versus the reference and there is no soft clipping.
#'
#' @return A a list with the elements:
#'
#' * \code{good.reads}, \code{data.frame} with
#' with column names for the
#' first 11 columns as specified in
#' https://en.wikipedia.org/wiki/SAM_(file_format)
#' with one row per read. The result contains only the
#' reads that meet certain conditions.
#'
#' * \code{high.flag.reads}, he same information as
#' \code{good.reads}, but for reads wth high FLAGs.
#'
#' * \code{low.mapq.reads}, the same information as
#' \code{good.reads}, but for reads wth low map quality.
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

  if (Sys.getenv("TRACEDBS") != "") {
    message("\nReadSamfile: filename = ", filename)
    message("ReadSamfile: starting with ", nrow(df), " reads")
  }

  df$FLAG <- suppressWarnings(as.numeric(df$FLAG))
  # remove weird rows if required
  bad.flag.reads <- df[is.na(df$FLAG)]
  df <- df[!is.na(df$FLAG), ]
  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: ",
            nrow(bad.flag.reads),
            " reads removed because FLAG could not interpreted as a number")
  }

  # Include only reads with FLAG < 256
  # Higher values mark reads that
  # - are "not pirmary alignment"
  # - failed vendor QC
  # - are PCR or optical duplicates
  # - supplementary alignments (e.g. split, split /inverted read)
  # https://broadinstitute.github.io/picard/explain-flags.html
  #
  bad.flag <- df$FLAG >= 256
  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: ",
            sum(bad.flag),
            " reads with bad FLAGs removed")
    if (sum(bad.flag) > 0) {
      message("bad FLAGs: ",
              paste(df[bad.flag, "FLAG"], collapse = ", "))
    }
  }
  df <- df[!bad.flag, ]
  bad.flag.reads <- rbind(bad.flag.reads, df[bad.flag, ])

  # Keep only reads that have a pair on the same chromosome.
  bad.mate.reads <- df[df$Mate_CHROM != "=", ]
  df <- df[df$Mate_CHROM == "=", ]
  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: ",
            nrow(bad.mate.reads),
            " reads removed pair was on another chromosome")
  }

  bad.mapq.reads <- df[df$MAPQ < 30, ]
  df <- df[df$MAPQ >= 30, ]
  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: ",
            nrow(bad.mapq.reads),
            " reads removed because MAPQ < 30")
  }

    if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: returning ", nrow(df), " reads")
  }

  return(df)
}
