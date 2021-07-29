PCAWG_prep_VCFs_and_BEDs_from_table <- function(table.name,
                                          out.dir,
                                          indiv.vcf.dir = "~/collab.vcf/",
                                          pcawg.vcf.dir = "~/pcawg.vcf/final_consensus_12aug_passonly/snv_mnv/"
                                          ) {
  tt <- data.table::fread(table.name)
  message("prep_PCAWG_VCFs_and_BEDs_from_table: ", nrow(tt), " rows to process")
  for (aliquot.id in tt$aliquot_id) {
    PCAWG_merge_callers(
      aliquot.id    = aliquot.id,
      indiv.vcf.dir = indiv.vcf.dir,
      pcawg.vcf.dir = pcawg.vcf.dir,
      out.dir       = out.dir)
  }
}

if (FALSE) {
  out.vcf.dir <- "~/mvv/short_test5"
  tt <-     "~/DBSverify/data-raw/short_collaboratory_bams.csv"
  PCAWG_prep_VCFs_and_BEDs_from_table(
    tt,
    out.dir = out.vcf.dir,
    indiv.vcf.dir = "~/collab.vcf/",
    pcawg.vcf.dir = "~/pcawg.vcf/final_consensus_12aug_passonly/snv_mnv/")

  # Then process the files after they have been downloaded from the collaboratory.
  # Sample code is in PCAWG_read_table_and_evaluate_DBS.R


}
