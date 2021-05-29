test_that("Long test of input/HepG2_AA1_20uM_SL_cl1_SNVresult.vcf -- full set of DBSs")
vcf.name <- "input/HepG2_AA1_20uM_SL_cl1_SNVresult.vcf"

vcf.list <-ICAMS::ReadAndSplitVCFs(vcf.name,
                                   variant.caller = "strelka",
                                   num.of.cores = 10,
                                   max.vaf.diff = 1)

DBS.vcf <- vcf.list$DBS[[1]]

evaluated.vcf <- VerifyDBSVcf(DBS.vcf,
                              Nbam = "input/HepG2_AA1_DBSlocs_Normal.bam",
                              Tbam = "input/HepG2_AA1_DBSlocs_Tumor.bam")

outfile = "xxx2.vcf"
write.table(evaluated.vcf, file = outfile, sep="\t", quote=F, row.names = F)
new <- data.table::fread(outfile)
old <- data.table::fread("input/HepG2_AA1_DBS_evaluated.vcf")
all.equal(new, old)
