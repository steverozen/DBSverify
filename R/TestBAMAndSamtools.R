TestBAMAndSamtools <- function(BAM.name, chrom.1.name = 1) {
  save.file.path <- tempfile()
  arb.pos <- 74823446
  CreateBAMSliceFileSamtools(BAM.name, chrom.1.name, arb.pos, padding = 10, save.file.path)
  sam <- ReadSamfile(save.file.path)
  if (nrow(sam) == 0) warning("No reads in ", BAM.name, " at ", chrom.1.name, ":", arb.pos)
  unlink(save.file.path)
  return(sam)
}
