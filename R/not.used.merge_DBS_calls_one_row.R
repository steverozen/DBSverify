not.used.merge_DBS_calls_one_row <- function(row, suffix.list) {
  row <- unlist(row)
  REF.cols <- paste("REF", suffix.list, sep = "")
  ALT.cols <- paste("ALT", suffix.list, sep = "")
  ref <- row[REF.cols]
  alt <- row[ALT.cols]
  ref.ok <- which(!is.na(ref))
  alt.ok <- which(!is.na(alt))
  ref2 <- unique(ref[ref.ok])
  alt2 <- unique(alt[alt.ok])
  stopifnot(length(ref2) == 1)
  stopifnot(length(alt2) == 1)
  rr <- c(row[c("CHROM", "POS")], REF = ref2, ALT = alt2, num.support = length(ref.ok))
  return(rr)
}
