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
  script.name <- paste0(vcf.name, "_igv_script.txt")

  devtools::load_all(".")
  debug(create_IGV_snapshot_script)
  create_IGV_snapshot_script(vcf.name = vcf.name,
                             Nbam.name = Nbam.name,
                             Tbam.name = Tbam.name,
                             out.dir   = file.path(getwd(),
                                                   "tests/testthat/snapshot"),
                             igv.script.name = script.name,
                             genome    = "Human hg38") # Different in Linux versus Windows


  igv.path <- "/home/gmssgr/igv/IGV_Linux_2.10.2/igv.sh"
  Sys.setenv(DISPLAY = "localhost:10.0")
  foo <- system2(igv.path, c(igv.path, "-b", script.name), env = c(DISPLAY = "localhost:10.0"))

  master.file <- "~/DBSverify/data-raw/short_collaboratory_bams.csv"

  tt <- data.table::fread(master.file)

  igv.path <- "/home/gmssgr/igv/IGV_Linux_2.10.2/igv.sh"
  debug(process.one.row)
  debug(create_IGV_snapshot_script)
  apply(X = tt, MARGIN = 1, process.one.row, igv.path = igv.path)

  }


process.one.row <- function(row, igv.path) {
  vcf.name <- paste0("/home/gmssgr/mvv/short_test5/",
                     row["icgc_donor_id"], "_", row["aliquot_id"], "_PCAWG_evaluated.vcf")
  stopifnot(file.exists(vcf.name))

  Tbam.name <- paste0("/home/gmssgr/mvv/bamSlice_folder/",
                      row["icgc_donor_id"], "_", row["T_Specimen ID"], "_dbs_srt_fromBED.bam")
  stopifnot(file.exists(Tbam.name))

  Nbam.name <- paste0("/home/gmssgr/mvv/bamSlice_folder/",
                      row["icgc_donor_id"], "_", row["N_Specimen ID"], "_dbs_srt_fromBED.bam")
  stopifnot(file.exists(Nbam.name))

  new.dir <- file.path("/home/gmssgr/mvv/short_test5/", row["icgc_donor_id"])

  igv.script.name <- file.path(new.dir, "igv_script.txt")

  create_IGV_snapshot_script(vcf.name = vcf.name,
                             Nbam.name = Nbam.name,
                             Tbam.name = Tbam.name,
                             out.dir   = new.dir,
                             igv.script.name = igv.script.name,
                             genome    = "Human hg19") # Different in Linux versus Windows

  foo <- system2(igv.path, args = c(igv.path, "-b", igv.script.name), wait = TRUE) # , env = c(DISPLAY = "localhost:10.0"))
  # convertGraph::convertGraph()
  # TraMineRextras::convert.g(path = new.dir, from = "pdf", to = "png")
  # qpdf::pdf_combine(input = dir(new.dir, pattern = "\\.*pdf"), output = path(new.dir, "combined.pdf"))
}


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
