# Basic test script; cannot use in production because cannot share BAMs

setwd("../mvv")

vcf.name <- "22.vcf"
Nbam <- "MELA-0161N_ch22.bam"
Tbam <- "MELA-0161T_ch22.bam"

xx <- ReadVCFAndBAMsAndProcess(
  vcf.name, Nbam.name = Nbam, Tbam.name = Tbam, variant.caller = "PCAWG")

old <- data.table::fread("MELA-0161_PCAWG_DBS_evaluated.vcf")
new <- data.table::fread(xx$vcf.name)

cat("regression result", all.equal(old, new), "\n")
