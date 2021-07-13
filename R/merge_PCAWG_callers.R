
merge_PCAWG_callers <- function(aliquot.id,
                                indiv.vcf.dir,
                                pcawg.vcf.dir,
                                out.dir = ".",
                                verbose = TRUE) {
  joined.dbs <-  join_PCAWG_callers(
    aliquot.id    = aliquot.id,
    indiv.vcf.dir = indiv.vcf.dir,
    pcawg.vcf.dir = pcawg.vcf.dir,
    verbose       = verbose)

  merged.dbs <-
    t(apply(
      joined.dbs,
      MARGIN = 1,
      FUN = merge_PCAWG_DBS_calls_one_row))

  merged.dbs <- data.table::data.table(merged.dbs)
  merged.dbs$POS <- as.numeric(merged.dbs$POS)
  colnames(merged.dbs) <-
    c("#CHROM", "POS", "REF", "ALT", "num.callers", "which.callers")

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
