# Test the case where different callers have different alternative alleles.
# This is a Collaboratory VCF, not sure if it is considered protected data,
# so the VCF is not in the repository or the R package.
PCAWG_prep_all_Collaboratory_VCFs_and_BEDs <- function() {
  out.vcf.dir <- "~/mvv/PCAWG_all_Collaboratory_VCFs_and_BEDs"
  tt          <-  "~/DBSverify/data-raw/13446fc4-9fa6-4bb3-a00a-b9265b41eb12.test.row.csv"
  prep_PCAWG_VCFs_and_BEDs_from_table(
    tt,
    out.dir = ".",
    indiv.vcf.dir = "~/collab.vcf/",
    pcawg.vcf.dir = "~/pcawg.vcf/final_consensus_12aug_passonly/snv_mnv/")
}
