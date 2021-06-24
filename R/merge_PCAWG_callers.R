
merge_PCAWG_callers <- function(aliquot.id, indiv.vcf.dir, pcawg.vcf.dir) {
  joined.dbs <-  join_PCAWG_callers(
    aliquot.id    = a.id1,
    indiv.vcf.dir = indiv.vcf.dir,
    pcawg.vcf.dir = pcawg.vcf.dir)

  merged.dbs <-
    t(apply(
      joined.dbs,
      MARGIN = 1,
      FUN = merge_DBS_calls_one_row,
      suffix.list = c("_mt", "_dk", "_ms", "_sa", "")))

  merged.dbs <- data.table::data.table(merged.dbs)
  merged.dbs$POS <- as.numeric(merged.dbs$POS)
  colnames(merged.dbs) <-
    c("#CHROM", "POS", "REF", "ALT", "num.callers")

  # This function writes the bed file
  bed <- VCF_to_BED(merged.dbs,
                    out.bed = paste0(aliquot.id, "_merged_PCAWG_DBS.bed"))

  data.table::fwrite(merged.dbs,
                     paste0(aliquot.id, "_merged_PCAWG_DBS.vcf"),
                     sep = "\t")

  return(invisible(merged.dbs))
}
