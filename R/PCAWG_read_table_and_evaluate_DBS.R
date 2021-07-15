#' Read a table that specifies how to process VCFs and mini-BAMS to evaluate DBS calls, for PCAWG Collaboratory data.
#'
#' @param in.table the file path of the table to process; in production,
#'   .../DBSverify/data-raw/collaboratory_bams_full_2021_07_13.csv,
#'   for testing .../DBSverify/data-raw/short_collaboratory_bams.csv.
#'
#' @param in.vcf.dir The path to the directory containing the DBS VCF files.
#'
#' @param minibam.dir The path to the directory containing the mini BAMs.
#'
#' @param out.vcf.dir The path to the directory in which to put the "evaluated"
#'   DBS VCF files.
#'
#' @details This is a specialized function for processing
#' PCAWG data from the "Collaboratory" once the
#' miniBAMs have been created in the Collaboratory and
#' downloaded. The \code{in.table} and associated BED files
#' were used to specify the contents of the miniBAMs.
#' The result consists of the "evaluated" DBS VCF files.
#' The naming of the input and output VCF files and the mini BAMs is governed
#' by the contents of \code{in.table}, with the VCF file names
#' incorporating the "aliquot_id" and the miniBAM names
#' based on the icgc_donor_id and the
#' T_Specimen ID and N_Specimen ID.
#'
#' @export

PCAWG_read_table_and_evaluate_DBS <- function(in.table,
                                              in.vcf.dir,
                                              minibam.dir,
                                              out.vcf.dir = in.vcf.dir) {

  tt <- data.table::fread(in.table)

  make.bam.name <- function(rrr, which.spec.id) {
    root <- paste0(rrr["icgc_donor_id"], "_", rrr[which.spec.id], "_dbs_srt_fromBED.bam")
    rx <- file.path(minibam.dir, root)
    return(rx)
  }

  Tbam.name <- apply(tt, MARGIN = 1, make.bam.name, "T_Specimen ID")

  Nbam.name <- apply(tt, MARGIN = 1, make.bam.name, "N_Specimen ID")


  in.vcf.name   <- lapply(tt$aliquot_id,
                          function(ai) {
                            return(
                              file.path(
                                in.vcf.dir,
                                paste0(ai, "_merged_PCAWG_DBS.vcf"))
                            )
                          })

  out.vcf.name   <- lapply(tt$aliquot_id,
                           function(ai) {
                             return(
                               file.path(
                                 in.vcf.dir,
                                 paste0(ai, "_PCAWG_evaluated.vcf"))
                             )
                           })

  for (ii in 1:nrow(tt)) {
    Read_DBS_VCF_and_BAMs_to_verify_DBSs(
      input.vcf     = in.vcf.name[[ii]],
      Nbam.name     = Nbam.name[[ii]],
      Tbam.name     = Tbam.name[[ii]],
      outfile       = out.vcf.name[[ii]],
      filter.status = NULL)
  }
}

if (FALSE) {


  devtools::load_all("~/DBSverify")
  debug(PCAWG_read_table_and_evaluate_DBS)
  # debug(Read_DBS_VCF_and_BAMs_to_verify_DBSs)
  # debug(ICAMS:::ReadVCFs)
  PCAWG_read_table_and_evaluate_DBS(in.table = "~/DBSverify/data-raw/short_collaboratory_bams.csv",
                                    in.vcf.dir  = "~/mvv/short_test5/",
                                    bamslice.dir = "~/mvv/bamSlice_folder/")




}
