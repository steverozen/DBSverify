ReadVCFAndBAMsAndProcess <- function(vcf.name, Nbam.name, Tbam.name, variant.caller) {

  # Add check for .bam.bai files here and in VerifyDBSVcf()

  # vcf.name <- "input/HepG2_AA1_20uM_SL_cl1_SNVresult.vcf"

  vcf.list <-ICAMS::ReadAndSplitVCFs(vcf.name,
                                     variant.caller = variant.caller,
                                     num.of.cores   = 10,
                                     max.vaf.diff   = 1)

  DBS.vcf <- vcf.list$DBS[[1]]

  evaluated.vcf <- VerifyDBSVcf(DBS.vcf, Nbam = Nbam.name, Tbam = Tbam.name)

  outfile <- paste0("evaluated_", vcf.name) # Maybe move the evaluated to the end of the file
  write.table(evaluated.vcf, file = outfile, sep="\t", quote=F, row.names = F)
}
