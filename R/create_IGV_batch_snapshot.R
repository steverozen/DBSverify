# Test script for creating IGV snapshots; creates a snapshot for each variant
# in each tumour.

# WARNING -- I cannot figure out how to specify the genome using the genome batch command.
#

if (FALSE) {

  master.file <- "~/DBSverify/data-raw/short_collaboratory_bams.csv"

  tt <- data.table::fread(master.file)

  igv.path <- "/home/gmssgr/igv/IGV_Linux_2.10.2/igv.sh"
  # debug(process.one.row)
  # debug(create_IGV_snapshot_script)
  apply(X = tt, MARGIN = 1, for.testing.process.one.row, igv.path = igv.path)

  }


for.testing.process.one.row <- function(row, igv.path) {
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
                             genome    = "hg19.genome") # Different in Linux versus Windows

  foo <- system2(igv.path, args = c(igv.path, "-b", igv.script.name), wait = TRUE) # , env = c(DISPLAY = "localhost:10.0"))
  # convertGraph::convertGraph()
  # TraMineRextras::convert.g(path = new.dir, from = "pdf", to = "png")
  # qpdf::pdf_combine(input = dir(new.dir, pattern = "\\.*pdf"), output = path(new.dir, "combined.pdf"))
}
