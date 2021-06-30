# get DBS and try to merge adjacent SBSs.  This has only been tested on the Hartwig data.

prep_DBS_and_adjacent_SBS <-
  function(vcf.name,
           remove.from.dbs = c("QUAL", "INFO", "FORMAT", "VAF", "read.depth"),
           remove.from.sbs = remove.from.dbs) {
    # vcfs <- ICAMS::ReadAndSplitMutectVCFs(vcf.name)
    vcfs <- ICAMS::ReadAndSplitVCFs(vcf.name,
                                    caller = "unknown")
    dbs.vcf <- vcfs$DBS[[1]]
    dbs.vcf <- remove_named_column(dbs.vcf, remove.from.dbs)

    sbs.vcf <- vcfs$SBS[[1]]

    sbs.vcf.file <- tempfile()
    colnames(sbs.vcf)[1] <- "#CHROM"
    sbs.vcf <- remove_named_column(sbs.vcf, remove.from.sbs)

    data.table::fwrite(data.table::as.data.table(sbs.vcf),
                       sbs.vcf.file,
                       sep       = "\t",
                       col.names = TRUE,
                       scipen    = 10000000) # Don't use scientific notation

    new.split <- ICAMS::ReadAndSplitVCFs(files           = sbs.vcf.file,
                                         variant.caller  = "unknown",
                                         always.merge.SBS = TRUE)

    new.dbs <- new.split$DBS[[1]]
    new.dbs <- remove_named_column(new.dbs,
                                   c("VAF", "read.depth", "remark.for.DBS"))

    colnames(new.dbs)[3:ncol(new.dbs)] <-
      paste(colnames(new.dbs)[3:ncol(new.dbs)], "ad", sep = "_")
    new.dbs$source.was.SBS <- TRUE

    j1 <- dplyr::full_join(dbs.vcf, new.dbs)

    j1[is.na(j1$source.was.SBS), "source.was.SBS"] <- FALSE

  return(invisible(j1))


}

if (FALSE) {
  prep_DBS_and_adjacent_SBS(
    "~/mvv/UTUC-test-case/AA_UTUC_Taiwan_T11_Mutect2_somatic_PASS.vcf")

  prep_DBS_and_adjacent_SBS("~/mvv/HMF_test.vcf")

  debug(prep_DBS_and_adjacent_SBS)
  # 717 DBS in original inpout
  # 8 more from adjacent SBSs
  foo <- prep_DBS_and_adjacent_SBS("~/mvv/HMF-VCFs/CPCT02020306T.purple.somatic.vcf.gz")


  }
