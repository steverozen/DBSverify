

ReadDBSVCF <- function(DBSVCF.name) {
  vcf <- utils::read.table(
    DBSVCF.name,
    header = TRUE,
    sep = "\t",
    colClasses = "character")
  vcf[ , "POS"] <- numeric(vcf[ , "POS"])
  return(vcf)
}
