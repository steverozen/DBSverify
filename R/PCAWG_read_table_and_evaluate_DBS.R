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
#' @param bam.suffix String to add to end of BAM file name; depends on
#'   the conventions used by the script (run on the Collaboratory)
#'   that generated the miniBAMs.
#'
#' @param verbose If > 0 generate some progress messages.
#'
#' @details This is a specialized function for processing
#' PCAWG data from the ICGC (International Cancer
#' Genome) "Collaboratory" cloud computing system, once the
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
                                              out.vcf.dir = in.vcf.dir,
                                              bam.suffix = "_dbs_srt",
                                              verbose     = 1) {

  tt <- data.table::fread(in.table)

  make.bam.name <- function(rrr, which.spec.id) {
    root <- paste0(rrr["icgc_donor_id"], "_", rrr[which.spec.id], bam.suffix, ".bam")
    rx <- file.path(minibam.dir, root)
    return(rx)
  }

  # Get a vector of the tumor BAM names based on tt.
  Tbam.name <- apply(tt, MARGIN = 1, make.bam.name, "T_Specimen ID")

  # Get a vector of the normal BAM names based on tt.
  Nbam.name <- apply(tt, MARGIN = 1, make.bam.name, "N_Specimen ID")

  # Get a vector of the input VCF file names based tt.
  in.vcf.name   <- lapply(tt$aliquot_id,
                          function(ai) {
                            return(
                              file.path(
                                in.vcf.dir,
                                paste0(ai, "_merged_PCAWG_DBS.vcf"))
                            )
                          })
  # Get a vector of the output VCF file names based on tt.
  out.vcf.name <-
    apply(tt,
          MARGIN = 1,
          FUN = function(rrr) {
            return(
              file.path(
                out.vcf.dir,
                paste0(rrr["icgc_donor_id"], "_", rrr["aliquot_id"],
                       "_PCAWG_evaluated.vcf"))
            )
          })

  for (ii in 1:nrow(tt)) {
    if (!CheckBAM(Nbam.name[[ii]], must.succeed = FALSE)) {
      message("Skipping ", in.vcf.name[[ii]], "; no corresponding normal BAM")
      next
    }
    if (!CheckBAM(Tbam.name[[ii]], must.succeed = FALSE)) {
      message("Skipping ", in.vcf.name[[ii]], "; no corresponding tumor BAM")
      next
    }
    try(
      Read_DBS_VCF_and_BAMs_to_verify_DBSs(
        input.vcf     = in.vcf.name[[ii]],
        Nbam.name     = Nbam.name[[ii]],
        Tbam.name     = Tbam.name[[ii]],
        outfile       = out.vcf.name[[ii]],
        filter.status = NULL,
        verbose       = verbose)
    )
  }
}

if (FALSE) {
  # For testing

  devtools::load_all("~/DBSverify")

  PCAWG_read_table_and_evaluate_DBS(in.table = "~/DBSverify/data-raw/short_collaboratory_bams.csv",
                                    in.vcf.dir  = "~/mvv/short_test5/",
                                    minibam.dir = "~/mvv/bamSlice_folder/",
                                    bam.suffix  = "_dbs_srt_fromBED")


  # Examining results

  ff <- function(filename) {
    nn <- data.table::fread(filename)
    print(unique(nn$DBSconclusion))
    View(nn)
  }

  ff("~/mvv/short_test5/DO52611_f82d213f-9ba5-7b6b-e040-11ac0c486882_PCAWG_evaluated.vcf")
  ff("~/mvv/short_test5/DO52610_f221c897-6ad0-0df9-e040-11ac0c4813ef_PCAWG_evaluated.vcf")
  ff("~/mvv/short_test5/DO52608_f8696c79-b165-92a6-e040-11ac0c4804bf_PCAWG_evaluated.vcf")
  ff("~/mvv/short_test5/DO52605_f82d213f-bc99-5b1d-e040-11ac0c486880_PCAWG_evaluated.vcf")
  ff("~/mvv/short_test5/DO52606_f856fa85-fdb8-c0b0-e040-11ac0d480b4e_PCAWG_evaluated.vcf")
  ff("~/mvv/short_test5/DO52612_f8467ec8-2d61-ba21-e040-11ac0c483584_PCAWG_evaluated.vcf")



}
