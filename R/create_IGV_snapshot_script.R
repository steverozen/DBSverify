create_IGV_snapshot_script <-
  function(vcf.name,
           Tbam.name,
           Nbam.name,
           igv.script.name = paste0(vcf.name, "_igv_script.txt"),
           out.dir,
           genome = "Human hg38",
           num.base.padding = 70) {

    # Do this early in case there is a problem with the VCF file
    vcf <- ICAMS::ReadVCFs(vcf.name, filter.status = NULL)[[1]]

    if (!dir.exists(out.dir)) {
      if (!dir.create(out.dir)) {
        warning("Did not create", out.dir)
      }
    }

    pp <- function(...) {
      cat(..., "\n", sep = "", file = igv.script.name, append = TRUE)
    }

    cat("new\n", file = igv.script.name)
    pp("snapshotDirectory ", out.dir)
    pp('genome "', genome, '"')
    pp("load ", '"', Tbam.name, '"')
    pp("load ", '"', Nbam.name, '"')


    for (ll in 1:nrow(vcf)) {
      ln <- vcf[ll, ]
      callers <- ln[["which.callers"]]
      if (!is.null(callers)) {
        callers <- gsub("[^A-Za-z0-9]", "_", callers)
      } else {
        callers <- ""
      }
      conc <- ln[["DBSconclusion"]]
      if (!is.null(conc)) {
        conc <- gsub("[^A-Za-z0-9]", "_", conc)
      } else {
        conc <- ""
      }
      pos <- ln[["POS"]]
      pp("goto ", ln[["CHROM"]], ":", pos - num.base.padding, "-", pos + num.base.padding)
      pp("sort base")
      pp("collapse")
      # pp("squish")
      pp("snapshot ", ln[["CHROM"]], "_", ln[["POS"]], "_", callers, "_", conc, ".png")
    }

    pp("exit")

  }
