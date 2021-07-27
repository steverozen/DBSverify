# Function to test creating IGV snapshots from cell line BAM files
# (so not protected data.)

test_create_IGV_snapshot_script <- function(
  igv.path = "/home/gmssgr/igv/IGV_Linux_2.10.2/igv.sh",
  run.IGV = FALSE,
  vcf.name

) {
  Nbam.name <- "tests/testthat/input/HepG2_AA1_DBSlocs_Normal.bam"
  Nbam.name <- file.path(getwd(), Nbam.name)
  Tbam.name <- "tests/testthat/input/HepG2_AA1_DBSlocs_Tumor.bam"
  Tbam.name <- file.path(getwd(), Tbam.name)
  script.name <- tempfile()
  out.dir <- tempfile()
  dir.create(out.dir)

  # devtools::load_all(".")
  # debug(create_IGV_snapshot_script)
  create_IGV_snapshot_script(vcf.name        = vcf.name,
                             Nbam.name       = Nbam.name,
                             Tbam.name       = Tbam.name,
                             out.dir         = out.dir,
                             igv.script.name = script.name,
                             # genome = "https://s3.amazonaws.com/igv.org.genomes/hg38/hg38.genome")
                             genome          = "Human (hg38)")
                             # genome          = "/home/gmssgr/igv/genomes/hg38.genome")

  if (run.IGV) {
    # Sys.setenv(DISPLAY = DISPLAY)
    foo <- system2(igv.path, args = c(igv.path, "-b", script.name), wait = TRUE)
  }
  rr <- scan(file = script.name, what = "character")
  return(rr)
}
