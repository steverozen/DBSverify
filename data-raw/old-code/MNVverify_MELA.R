# Basic test scritp

setwd("../mvv")

smp <- "MELA-0161"
vcf.name <- "22.vcf"
Nbam <- "MELA-0161N_ch22.bam"
Tbam <- "MELA-0161T_ch22.bam"


vcf.list <-ICAMS::ReadAndSplitVCFs(vcf.name,
                 variant.caller = "unknown",
                 num.of.cores = 10,
                 names.of.VCFs = smp,
                 get.vaf.function = GetPCAWGVAF,
                 max.vaf.diff = 1)

DBS.vcf <- vcf.list$DBS[[1]]

evaluated.vcf <- VerifyDBSVcf(DBS.vcf, Nbam = Nbam, Tbam = Tbam)

outfile <- paste0(smp,"_PCAWG_DBS_evaluated.vcf")
write.table(evaluated.vcf, file = outfile, sep="\t", quote=F, row.names = F)

old <- data.table::fread("MELA-0161_PCAWG_DBS_evaluated.vcf")
new <- data.table::fread(outfile)
all.equal(old, new)
cat("regression result", all.equal(old, new), "\n")

