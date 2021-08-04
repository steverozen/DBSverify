#' Read a "SAM file", discarding some reads that cannot be interpreted for our purposes.
#'
#' @details
#' SAM stands for "Sequence Alignment Map", a text file that
#' represents aligned next-generation sequencing reads
#' (https://en.wikipedia.org/wiki/SAM_(file_format)).
#' Only keep reads that satisfy certain conditions:
#'
#' * Mate pair maps to same chromosome
#' * Mapping quality >= 30
#' * "FLAG" < 256 (see info on the return value
#' element \code{reads.with.bad.FLAG}, below).
#' * The CIGAR string
#'  is only \code{\\d+M} (one or more digits followed M, and nothing
#'  else before or after).
#'  This means there are no insertions or deletions in the read
#'  versus the reference and there is no soft clipping. Insertions
#'  or deletions or clipping shift the read within the
#'  SAM file, after which this function cannot keep track of the
#'  DBS location.
#'
#' @param filename The name of the SAM file to read.
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
#' * \code{total.depth}, the initial depth including
#' "bad" reads
#'
#' * \code{reads.with.bad.FLAG}. These are reads
#'  with FLAG >= 256, which marks reads that
#' (i) are "not primary alignment"
#' (ii) failed vendor QC
#' (iii) are PCR or optical duplicates
#' (iv) are supplementary alignments (e.g. split, split /inverted read).
#'  See https://broadinstitute.github.io/picard/explain-flags.html
#'
#' * \code{reads.with.bad.CIGAR}, reads with CIGAR string
#'  that indicates an indel in the read or clipping.
#'
#' * \code{reads.with.bad.MAPQ}, reads with MAPQ < 30.
#'
#' * \code{reads.with.bad.Mate_CHROM}, reads with a
#'   mate on a different chromosome.

ReadSamfile <- function(filename) {
  sam.lines <- readLines(filename)

  # SAM file headers start with "@"; discard them and keep the
  # lines that contains reads.
  read.lines <- sam.lines[-grep("^@", sam.lines)]

  # Perfectly duplicated reads arise during generation of miniBAMs when
  # reads overlap more than one DBS position. We need to
  # remove these.
  dups <- duplicated(read.lines)
  if (any(dups)) {
    # For development we log this information.
    message("dup reads in SAM file ", filename)
    read.lines <- read.lines[!dups]
  }

  possible.good.reads <-
    as.data.frame(stringr::str_split_fixed(read.lines, "\t", 16))
  colnames(possible.good.reads)[1:11] <-
    c("QNAME", "FLAG", "CHROM", "POS", "MAPQ", "CIGAR", "Mate_CHROM",
      "Mate_POS", "InsertSize", "SEQ", "QUAL")

  if (Sys.getenv("TRACEDBS") != "") {
    message("\nReadSamfile: filename = ", filename)
    message("ReadSamfile: starting with ", nrow(possible.good.reads), " reads")
  }

  possible.good.reads$FLAG <-
    suppressWarnings(as.numeric(possible.good.reads$FLAG))
  # remove weird rows if required
  bad.flags <- is.na(possible.good.reads$FLAG)
  reads.with.bad.FLAG <- possible.good.reads[bad.flags, ]
  possble.good.reads <- possible.good.reads[!bad.flags, ]
  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: ",
            nrow(reads.with.bad.FLAG),
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
  bad.flags <- possible.good.reads$FLAG >= 256
  reads.with.bad.FLAG <- rbind(reads.with.bad.FLAG, possible.good.reads[bad.flags, ])
  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: ",
            sum(bad.flags),
            " reads with bad FLAGs removed")
    if (sum(bad.flags) > 0) {
      message("bad FLAGs: ",
              paste(reads.with.bad.FLAG[ , "FLAG"], collapse = ", "))
    }
  }
  possible.good.reads <- possible.good.reads[!bad.flags, ]

  cigar.good <- grep("^\\d+M$", possible.good.reads$CIGAR)
  reads.with.bad.CIGAR <- possible.good.reads[-cigar.good, ]
  possible.good.reads <- possible.good.reads[cigar.good, ]
  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: ",
            nrow(reads.with.bad.CIGAR),
            " reads removed because CIGAR string was not \\d+M")
    message("Bad CIGARs: ", paste(reads.with.bad.CIGAR$CIGAR, collapse = ", "))
  }

  bad.mapq <- possible.good.reads$MAPQ < 30
  reads.with.bad.MAPQ <- possible.good.reads[bad.mapq, ]
  possible.good.reads <- possible.good.reads[!bad.mapq, ]
  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: ",
            nrow(reads.with.bad.MAPQ),
            " reads removed because MAPQ < 30")
  }

  # Keep only reads that have a pair on the same chromosome.
  ok.mate <- possible.good.reads$Mate_CHROM == "="
  reads.with.bad.Mate_CHROM <- possible.good.reads[!ok.mate, ]
  possible.good.reads       <- possible.good.reads[ok.mate, ]
  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: ",
            nrow(reads.with.bad.Mate_CHROM),
            " reads removed because mate was on another chromosome")
  }

  if (Sys.getenv("TRACEDBS") != "") {
    message("ReadSamfile: returning ", nrow(possible.good.reads), " reads")
  }

  rr <- list(good.reads = possible.good.reads,
             reads.with.bad.FLAG       = reads.with.bad.FLAG,
             reads.with.bad.CIGAR      = reads.with.bad.CIGAR,
             reads.with.bad.MAPQ       = reads.with.bad.MAPQ,
             reads.with.bad.Mate_CHROM = reads.with.bad.Mate_CHROM
             )

  # return(possible.good.reads)
  return(rr)

}
