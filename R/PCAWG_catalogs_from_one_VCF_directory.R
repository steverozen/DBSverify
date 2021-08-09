PCAWG_catalogs_from_one_VCF_directory <-
  function(vcf.dir, output.dir = vcf.dir) {
    old.warn <- getOption("warn")
    on.exit(options(warn = old.warn))
    options(warn = 1)
    vcfs <- dir(vcf.dir, full.names = TRUE, pattern = "\\.vcf")
    vcf.aliquot.id1 <- gsub(vcf.dir, "", vcfs)
    vcf.aliquot.id2  <- gsub("/DO\\d+_", "", vcf.aliquot.id1, perl = TRUE)
    vcf.aliquot.id  <- gsub("SP\\d+_", "", vcf.aliquot.id2, perl = TRUE)
    aliquot.id <- gsub("_PCAWG_evaluated.vcf", "", vcf.aliquot.id)
    tt <- data.table::fread(
      "~/DBSverify/data-raw/production_scripts/collaboratory_bams_2021_07_16.csv")

    jj <- dplyr::left_join(
      data.table::data.table(aliquot_id = aliquot.id), tt)

    xx <- lapply(vcfs, PCAWG_only_true_DBSs)

    # Sketch -- generate some kind of stats

    final.vcf.list <- lapply(xx, `[[`, "final.vcf")
    pcawg.vcf.list <- lapply(xx, `[[`, "pcawg.vcf")
    num.orig.dbs <- unlist(lapply(xx, `[[`, "orig.pcawg.dbs"))
    num.good.dbs <- unlist(lapply(xx, `[[`, "good.n.dbs"))

    my.names <- PCAWG7:::map_SP_ID_to_tumor_type(jj$`T_Specimen ID`) # Check that map_... exported,then fix :::
    names(final.vcf.list) <- my.names
    names(pcawg.vcf.list) <- my.names
    deltas <- (num.orig.dbs - num.good.dbs)
    names(deltas) <- my.names
    cc <- cbind(my.names, num.orig.dbs, num.good.dbs, deltas, num.good.dbs / num.orig.dbs)
    colnames(cc) <- c("PCAWG_ID", "Num_orig_PCAWG_DBS", "Num_good_DBS", "Difference", "Ratio")
    cc <- cc[order(as.numeric(cc[, "Difference"]), decreasing = TRUE), ]

    dbs.cats <- ICAMS::VCFsToDBSCatalogs(
      final.vcf.list,
      ref.genome = "hg19",
    )

    c78 <- dbs.cats$catDBS78
    c78 <- c78[ , order(colnames(c78))]

    ICAMS::WriteCatalog(c78,
                        file.path(output.dir, "DBS78.csv"))
    ICAMS::PlotCatalogToPdf(c78,
                            file.path(output.dir, "DBS78.pdf"))

    data.table::fwrite(data.table::data.table(cc),
                       file.path(output.dir, "summary.csv"))
  }

if (FALSE) {
  # ask nanhai to debug
all.output <- ICAMS::VCFsToCatalogs(
  files = vcfs,
  ref.genome = "hg19",
  num.of.cores = 1,
  region = "genome"
)
}

