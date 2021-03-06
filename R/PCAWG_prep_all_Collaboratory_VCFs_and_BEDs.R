# Run once to generate all VCFs and BEDs for the PCAWG Collabortory Data
# The file locations are specific to the locations when the data were generated
# The VCFs and BED files are protected data.
PCAWG_prep_all_Collaboratory_VCFs_and_BEDs <- function(tt) {
  out.vcf.dir <- "~/mvv/PCAWG_all_Collaboratory_VCFs_and_BEDs"
  PCAWG_prep_VCFs_and_BEDs_from_table(
    tt,
    out.dir = out.vcf.dir,
    indiv.vcf.dir = "~/collab.vcf/",
    pcawg.vcf.dir = "~/pcawg.vcf/final_consensus_12aug_passonly/snv_mnv/")
}
