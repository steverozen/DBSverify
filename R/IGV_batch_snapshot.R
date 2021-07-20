if (FALSE) {
  igv <- file.path("c:", "Program Files", "IGV_2.8.13", "igv.bat")

  shell("java --module-path=lib -Xmx4g @igv.args --module=org.igv/org.broad.igv.ui.Main")

  shell("java --module-path=lib -Xmx4g --module=org.igv/org.broad.igv.ui.Main")

  shell.exec(igv)

  bamfolder <-
    file.path("c:", "Users", "steve rozen", "Desktop", "bamSlice_folder")

  Tbam.name <-
    dir(file.path(bamfolder, "DO52605_SP116496_dbs_srt_fromBED.bam"),
        full.names = TRUE)

  Nbam.name <-
    dir(file.path(bamfolder, "DO52605_SP116498_dbs_srt_fromBED.bam"),
        full.names = TRUE)

  Tbam.name <-
    dir(file.path(bamfolder, "DO52605_SP116496_dbs_srt_fromBED.bam"),
        full.names = TRUE)

  Nbam.name <- "tests/testthat/input/HepG2_AA1_DBSlocs_Normal.bam"
  Nbam.name <- file.path(getwd(), Nbam.name)
  Tbam.name <- "tests/testthat/input/HepG2_AA1_DBSlocs_Tumor.bam"
  Tbam.name <- file.path(getwd(), Tbam.name)
  vcf.name <- "tests/testthat/input/HepG2_AA1_DBS_evaluated.vcf" # or HepG2_AA1.vcf.evaluated.vcf.regress?

  devtools::load_all(".")
  debug(create_IGV_snapshot_script)
  create_IGV_snapshot_script(vcf.name = vcf.name,
                             Nbam.name = Nbam.name,
                             Tbam.name = Tbam.name,
                             out.dir   = file.path(getwd(),
                                                   "tests/testthat/snapshot"))



}


create_IGV_snapshot_script <-
  function(vcf.name,
           Tbam.name,
           Nbam.name,
           igv.script.name = paste0(vcf.name, "_igv_script.txt"),
           out.dir,
           genome = "Human hg19",
           num.base.padding = 50) {

    # Do this early in case there is a problem with the VCF file
    vcf <- ICAMS:::ReadVCF(vcf.name)

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
    pp("genome ", genome)
    pp("load ", '"', Tbam.name, '"')
    pp("load ", '"', Nbam.name, '"')


    for (ll in 1:nrow(vcf)) {
      ln <- vcf[ll, ]
      pos <- ln[["POS"]]
      pp("goto ", ln[["CHROM"]], ":", pos - num.base.padding, "-", pos + num.base.padding)
      pp("sort base")
      pp("collapse")
      pp("snapshot ", ln[["CHROM"]], "_", ln[["POS"]], ".png")
    }

  }
