PCAWG_only_true_DBSs <- function(VCF.file) {
  vcf <- data.table::fread(VCF.file)
  colnames(vcf)[1] <- "CHROM"
  vcf$CHROM <- as.character(vcf$CHROM)
  orig.n.dbs <- nrow(vcf)
  pcawg.vcf <- vcf[grep("pcawg", vcf$which.callers), ]
  orig.pcawg.dbs <- nrow(pcawg.vcf)
  final.vcf <- pcawg.vcf[pcawg.vcf$DBSconclusion == "True DBS", ]
  return(list(final.vcf      = final.vcf,
              pcawg.vcf      = pcawg.vcf,
              orig.n.dbs     = orig.n.dbs,
              orig.pcawg.dbs = orig.pcawg.dbs,
              good.n.dbs     = nrow(final.vcf)))
}
