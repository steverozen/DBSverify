# Join info on BAMs from the ICGC Collaboratory and VCF
# PCAWG consensus VCF files (which are named by the aliquot
# ID of the tumor sample), so we know which
# VCF file goes with which BAM file names and BAM object ids.

# Run this function once to create the merged fused table.

PCAWG_fuse_sample_sheet_and_collab_manifest <- function() {

  # This is the "sample sheet" from PCAWG project.
  # Downloaded from
  #
  #   https://dcc.icgc.org/api/v1/download?fn=/PCAWG/data_releases/latest/pcawg_sample_sheet.v1.4.2016-09-14.tsv
  #
  # on 9 June 2021
  samp.sheet <- data.table::fread("data-raw/pcawg_sample_sheet.v1.4.2016-09-14.tsv")

  samp.sheet <- samp.sheet[samp.sheet$library_strategy == "WGS", ] # 5789 rows

  # This is the info from the Collaborary, using the query URL:
  #   https://dcc.icgc.org/repositories?filters=%7B%22file%22:%7B%22repoName%22:%7B%22is%22:%5B%22Collaboratory%20-%20Toronto%22%5D%7D,%22dataType%22:%7B%22is%22:%5B%22Aligned%20Reads%22%5D%7D,%22experimentalStrategy%22:%7B%22is%22:%5B%22WGS%22%5D%7D,%22study%22:%7B%22is%22:%5B%22PCAWG%22%5D%7D,%22fileFormat%22:%7B%22is%22:%5B%22BAM%22%5D%7D,%22software%22:%7B%22is%22:%5B%22BWA%20MEM%22,%22_missing%22%5D%7D%7D%7D&files=%7B%22from%22:1,%22size%22:25%7D
  #
  # The query URL gives us info on PCAWG BAM files stored at the Collaboratory.
  #
  # Then need to hover over the little organge page icon on the far right, just under
  # the text "Save/Edit File Results", which gives us the option "Export table as TSV".
  #
  # Downloaded 9 June 2021, 9:00 PM Singapore time

  bam.info <- data.table::fread(
    "data-raw/collaboratory_WGS_BAMS_repository_1623246282.tsv") # 7208 rows

  stopifnot(unique(bam.info$format) == "BAM")

  bam.info <- bam.info[  ,
                         c("Object ID",
                           "File Name",
                           "ICGC Donor",
                           "Specimen ID",
                           "Specimen Type",
                           "Size (bytes)")]

  si <- split(bam.info, f = bam.info$`ICGC Donor`) # 1830 elements (donors)

  organize.one.donor <- function(ssi) {
    # ssi is a table with ICGG Donor (e.g. DO9999) and Specimen ID (e.g. SP9999)
    # and indication of which are normals and which are tumor; there are full
    # BAMs and mini BAMS in separate rows.

    # For each Specimen ID get the row for the largest BAM (i.e. the full BAM).
    by.spec <- split(ssi, ssi$`Specimen ID`)
    max.specs <- lapply(by.spec, function(one.spec) {
      max.idx <- which(one.spec$`Size (bytes)` == max(one.spec$`Size (bytes)`))
      if (length(max.idx) > 1) stop("max not unique")
      return(one.spec[max.idx, ])
    })

    # max.specs is a list with one elment for each Specimen ID.
    return(max.specs)
  }

  donor.pairs <- lapply(si, organize.one.donor)
  # donor.pairs is a nested list structure with one top level element for each
  # donor, one sub element for the normal, and one or more sub elements for the
  # tumor samples.

  # For sanity check; these are donors with more than one tumor sample.
  print(which(unlist(lapply(donor.pairs, length)) > 2))
  # 302, 631, ...

  no.match <- which(unlist(lapply(donor.pairs, length)) == 1)
  print("One donor had no matched normal: ")
  print(no.match)
  donor.pairs <- donor.pairs[-no.match]

  is.normal <- function(one.specimen) {
    return(grepl("Normal", one.specimen$`Specimen Type`))
  }

  # Create one or more rows of tumor-normal pairs; all rows
  # have the same normal sample.
  flatten.one.donor <- function(pair) {
    norm.index <- which(unlist(lapply(pair, is.normal)))
    if (length(norm.index) == 0) {
      browser()
      stop("No normal sample")
    }
    pp <- pair[-norm.index]
    ppp <- lapply(pp, merge.tumor.and.normal, pair[[norm.index]])
    ppp <- do.call(rbind, ppp)
    return(ppp)
  }

  # Create a new row with info from one tumor and one normal.
  merge.tumor.and.normal <- function(tumor.row, normal.row) {
    colnames(tumor.row) <- paste0("T_", colnames(tumor.row))
    colnames(normal.row) <- paste0("N_", colnames(normal.row))
    rr <- cbind(tumor.row, normal.row)
    return(rr)
  }

  tumor.normal.pair <- tibble::as_tibble(do.call(rbind, lapply(donor.pairs, flatten.one.donor)))

  jj <- dplyr::inner_join(samp.sheet, tumor.normal.pair, by  = c("icgc_specimen_id" = "T_Specimen ID"))

  # Write a table with lots of extraneous information (probably useless).
  data.table::fwrite(jj, "data-raw/collaboratory_bams_full_2021_07_13.csv")

  ## Simplify the columns
  jjj <- jj[ , c("aliquot_id", "icgc_donor_id", "icgc_specimen_id",
                 "T_Object ID", "T_File Name",
                 "N_Specimen ID", "N_Object ID", "N_File Name")]
  colnames(jjj)[3] <- "T_Specimen ID"

  data.table::fwrite(jjj, "data-raw/collaboratory_bams_2021_07_13.csv")
}
