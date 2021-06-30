prep_PCAWG_VCFs_and_BEDs_from_table <- function(table.name,
                                          out.dir = "~/mvv/short_collab_test_2/",
                                          indiv.vcf.dir = "~/collab.vcf/",
                                          pcawg.vcf.dir = "~/pcawg.vcf/final_consensus_12aug_passonly/snv_mnv/"
                                          ) {
  tt <- data.table::fread(table.name)

  for (aliquot.id in tt$aliquot_id) {
    merge_PCAWG_callers(
      aliquot.id    = aliquot.id,
      indiv.vcf.dir = indiv.vcf.dir,
      pcawg.vcf.dir = pcawg.vcf.dir,
      out.dir       = out.dir)
  }
}

if (FALSE) {


  prep_PCAWG_VCFs_and_BEDs_from_table(
    "./data-raw/short_collaboratory_bams.csv",
    out.dir = "~/mvv/short_test3",
    indiv.vcf.dir = "~/collab.vcf/",
    pcawg.vcf.dir = "~/pcawg.vcf/final_consensus_12aug_passonly/snv_mnv/")


}
