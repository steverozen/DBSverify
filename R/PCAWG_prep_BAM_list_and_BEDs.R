#' Specialized function for processing PCAWG Collaboratory BAMS.
#'
#' This function is unfinished; it would need to be completed to generate
#' a table with bed files and instructions for miniBAMs for tumor BAMS.
#' However, this is not necessary as we already have the
#' necessary miniBAMs for donors with single tumor samples.
#'
#' @param only.multiple.tumors.per.normal If \code{TRUE} only generate output
#'   for tumors from donors with multiple tumor VCF? files.
#'
#' @keywords internal

PCAWG_prep_BAM_list_and_BEDs <- function(only.multiple.tumors.per.normal = TRUE) {
  tt <- data.table::fread("~/DBSverify/data-raw/production_scripts/collaboratory_bams_2021_07_16.csv")
  xx <- split(tt, tt$icgc_donor_id)
  input.vcf.dir <- "~/mvv/PCAWG_all_Collaboratory_VCFs_and_BEDs/"
  output.vcf.dir <- file.path(input.vcf.dir, "patch_dups")

  make.tumor.beds.and.table <- function(tt.row) {
    # Sketch only
    # read tumor vcf file
    # make tumor bed file
    # write row with tumor bed file name, tumor bam file id, tumor bam file name
  }

  if (FALSE) {
    lapply(tt, make.tumor.beds.and.table)
  }

  evaluate.info.table <- "DBSverify.info.csv"
  cat("aliquot.id,icgc.donor.id,T.specimen.id,N.specimen.id\n",
      file = evaluate.info.table)
  get.minibam.info.table <- "get.miniBAM.info.csv"
  cat("NBAM.bed.name,nbam.id,nbam.name\n",
      file = get.minibam.info.table)

  make.normal.beds.and.table <- function(xx.row) {
   if ((nrow(xx.row) == 1) && only.multiple.tumors.per.normal) return(NULL);
   # View(xx.row)
    # browser()
    n.spec.id <- unique(xx.row$`N_Specimen ID`)
    nbam.id   <- unique(xx.row$`N_Object ID`)
    nbam.name <- unique(xx.row$`N_File Name`)
    donor.id  <- unique(xx.row$icgc_donor_id)

    for (ii in 1:nrow(xx.row)) {
      aliquot.id <- xx.row[ii, "aliquot_id"]
      T.specimen.id <- xx.row[ii, "T_Specimen ID"]

      cat(paste(aliquot.id, donor.id, T.specimen.id, n.spec.id,
                sep = ","),
          "\n", sep = "", file = evaluate.info.table, append = TRUE)

      in.vcf.file <- file.path(input.vcf.dir,
                               paste0(aliquot.id, "_merged_PCAWG_DBS.vcf"))
      vcf <- data.table::fread(in.vcf.file)
      vcf$start   <- vcf$POS - 10
      vcf$end     <- vcf$POS + 11
      colnames(vcf)[1] <- "chrom"
      vcf.granges <- GenomicRanges::makeGRangesFromDataFrame(
        vcf,
        ignore.strand = TRUE)
      if (ii == 1) {
        final.granges <- vcf.granges
      } else {
        final.granges <-
          GenomicRanges::reduce(
            unlist(
              GenomicRanges::GRangesList(vcf.granges, final.granges)))
      }
    }
    NBAM.bed.name <- paste0(donor.id, "_", n.spec.id, "_NBAM.bed")
    df <- as.data.frame(final.granges)
    data.table::fwrite(
      df[ , c(1:3)],
      col.names = FALSE,
      sep = " ",
      file = file.path(output.vcf.dir, NBAM.bed.name)                                      )

    cat(paste(NBAM.bed.name, nbam.id, nbam.name, sep = ","),
        "\n", sep = "", file = get.minibam.info.table, append = TRUE)

  }

  lapply(xx, make.normal.beds.and.table)

}
