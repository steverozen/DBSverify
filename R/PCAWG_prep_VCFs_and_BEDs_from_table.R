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

  # Now process the files after they have been downloaded from the collaboratory
  bam.dir <- "~/mvv/collab_minibam"
  rows <- data.table::fread(tt)
  for (rrr in 1:nrow(rows)) {
    rr <- rows[rrr, ]
    N.bam.root  <- paste(rr$icgc_donor_id, rr$`N_Specimen ID`, "dbs.sorted.bam", sep = "_")
    T.bam.root  <- paste(rr$icgc_donor_id, rr$`T_Specimen ID`, "dbs.sorted.bam", sep = "_")
    DBSverify::Read_DBS_VCF_and_BAMs_to_verify_DBSs(
      input.vcf = file.path(out.vcf.dir, paste0(rr$aliquot_id, "_merged_PCAWG_DBS.vcf")),
      Nbam.name = file.path(bam.dir, N.bam.root),
      Tbam.name = file.path(bam.dir, T.bam.root),
      verbose   = 1,
      outfile   = file.path(out.vcf.dir, paste0(rr$aliquot_id, "_evaluated.vcf")),
      filter.status = NULL
    )
  }

  View(data.table::fread(file.path(out.vcf.dir, "f221c897-6ad0-0df9-e040-11ac0c4813ef_evaluated.vcf")))
  View(data.table::fread(file.path(out.vcf.dir, "f82d213f-bc99-5b1d-e040-11ac0c486880_evaluated.vcf")))


}
