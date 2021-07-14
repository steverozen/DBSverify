PCAWG_merge_DBS_calls_one_row <- function(row) {
  suffix.list <- c("_mt", "_dk", "_ms", "_sa", "_pcawg")
  row <- unlist(row)
  REF.cols <- paste("REF", suffix.list, sep = "")
  ALT.cols <- paste("ALT", suffix.list, sep = "")
  ref <- row[REF.cols]
  alt <- row[ALT.cols]
  ref.ok <- which(!is.na(ref))
  which.caller <- paste(suffix.list[ref.ok], collapse = ";")
  alt.ok <- which(!is.na(alt))
  ref2 <- unique(ref[ref.ok])
  alt2 <- unique(alt[alt.ok])
  if (length(ref2) != 1) {
    warning("Row with more than one reference allele: ",
            paste(ref2, collapse = ", "))
    bad.row <- utils::capture.output(print(row))
    warning(bad.row)
    rev2 <- row["REF_pcawg"]
  }
  if (length(alt2) != 1) {
    warning("Row with more than one alternate allele: ",
            paste(alt2, collapse = ", "))
    bad.row <- utils::capture.output(print(row))
    warning(bad.row)
    alt2 <- row["ALT_pcawg"]
  }
  rr <- c(row[c("CHROM", "POS")], REF = ref2, ALT = alt2, num.callers = length(ref.ok), which.caller = which.caller)
  return(rr)
}
