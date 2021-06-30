prep_DBS_and_adjacent_SBS <- function(vcf.name) {
  # vcfs <- ICAMS::ReadAndSplitMutectVCFs(vcf.name)
  vcfs <- ICAMS::ReadAndSplitVCFs(vcf.name,
                                  caller = "unknown")
  dbs.vcf <- vcfs$DBS[[1]]
  sbs.vcf <- vcfs$SBS[[1]]

  sbs.vcf.file <- tempfile()
  colnames(sbs.vcf)[1] <- "#CHROM"
  sbs.vcf <- remove_named_column(sbs.vcf, "VAF")

  data.table::fwrite(data.table::as.data.table(sbs.vcf),
                     sbs.vcf.file,
                     sep       = "\t",
                     col.names = TRUE,
                     scipen    = 10000000) # Don't use scientific notation

  new.split <- ICAMS::ReadAndSplitVCFs(files           = sbs.vcf.file,
                                       variant.caller  = "unknown",
                                       always.merge.SBS = TRUE)

  new.dbs <- new.split$DBS[[1]]
  new.dbs <- remove_named_column(new.dbs, "VAF")

  colnames(new.dbs)[3:ncol(new.dbs)] <-
    paste(colnames(new.dbs)[3:ncol(new.dbs)], "ad", sep = "_")

  j1 <- dplyr::full_join(dbs.vcf, new.dbs)

  return(0)


}

if (FALSE) {
  prep_DBS_and_adjacent_SBS(
    "~/mvv/UTUC-test-case/AA_UTUC_Taiwan_T11_Mutect2_somatic_PASS.vcf")

  prep_DBS_and_adjacent_SBS("~/mvv/HMF_test.vcf")

  debug(prep_DBS_and_adjacent_SBS)
  # 717 DBS in original inpout
  # 8 more from adjacent SBSs
  prep_DBS_and_adjacent_SBS("~/mvv/HMF-VCFs/CPCT02020306T.purple.somatic.vcf.gz")


  }
