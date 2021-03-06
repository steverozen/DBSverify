
PCAWG_merge_callers <- function(aliquot.id,
                                indiv.vcf.dir,
                                pcawg.vcf.dir,
                                out.dir = ".",
                                verbose = TRUE) {
  joined.dbs <-  PCAWG_join_callers(
    aliquot.id    = aliquot.id,
    indiv.vcf.dir = indiv.vcf.dir,
    pcawg.vcf.dir = pcawg.vcf.dir,
    verbose       = verbose)

  if (nrow(joined.dbs) == 0) {
    warning("No DBSs in ", aliquot.id)
    return(data.frame("#CHROM" = 0, POS = 0, ALT = 0, num.callers = 0, which.callers = 0)[-1, ])
  }

  merged.dbs <-
    t(apply(
      joined.dbs,
      MARGIN = 1,
      FUN = PCAWG_merge_DBS_calls_one_row))

  merged.dbs <- data.table::data.table(merged.dbs)
  colnames(merged.dbs) <-
    c("#CHROM", "POS", "REF", "ALT", "num.callers", "which.callers")
  merged.dbs$POS <- as.numeric(merged.dbs$POS)

  merged.dbs <- merged.dbs[!is.na(merged.dbs$ALT), ]
  merged.dbs <- merged.dbs[!is.na(merged.dbs$REF), ]

  if (verbose) message("Writing bed file")
  bed <- VCF_to_BED(
    merged.dbs,
    out.bed = file.path(out.dir, paste0(aliquot.id, "_merged_PCAWG_DBS.bed")))

  if (verbose) {
    message("Writing vcf file: ", nrow(merged.dbs), " rows")
  }


  data.table::fwrite(merged.dbs,
                     file.path(out.dir,
                               paste0(aliquot.id, "_merged_PCAWG_DBS.vcf")),
                     sep = "\t")

  return(invisible(merged.dbs))
}
