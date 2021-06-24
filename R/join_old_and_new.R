join_old_and_new_vcf <- function(old, new) {
  if (colnames(old)[1] == "#CHROM") {
    old.new <- dplyr::full_join(
      old,
      new,
      by = c("#CHROM" = "#CHROM", "POS" = "POS"))
    xold.new <- old.new[ ,

                         c("#CHROM", "POS", "NreadSupport.x", "TreadSupport.x", "NreadSupport.y", "TreadSupport.y", "DBSconclusion.x", "DBSconclusion.y")]

  } else if (colnames(old)[1] == "CHROM") {
    old.new <- dplyr::full_join(
      old,
      new,
      by = c("CHROM" = "CHROM", "POS" = "POS"))
    xold.new <-
      old.new[ , c("CHROM", "POS",
                   "NreadSupport.x", "TreadSupport.x",
                   "NreadSupport.y", "TreadSupport.y",
                   "DBSconclusion.x", "DBSconclusion.y")]
  } else {
    stop("First column name is not CHROM or #CHROM")
  }
  return(xold.new)
}
