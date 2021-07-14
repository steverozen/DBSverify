# Test the case where different callers have different alternative alleles.
# This is a Collaboratory VCF, not sure if it is considered protected data,
# so the VCF is not in the repository or the R package.
PCAWG_prep_all_Collaboratory_VCFs_and_BEDs_cd <- function() {
  out.vcf.dir <- "~/mvv/PCAWG_all_Collaboratory_VCFs_and_BEDs"
  tt          <-  "~/DBSverify/data-raw/4e7cdeda-6dc1-4f17-b853-72a68e5aa7e1.test.row.csv"
  prep_PCAWG_VCFs_and_BEDs_from_table(
    tt,
    out.dir = ".",
    indiv.vcf.dir = "~/collab.vcf/",
    pcawg.vcf.dir = "~/pcawg.vcf/final_consensus_12aug_passonly/snv_mnv/")
}
