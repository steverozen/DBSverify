# Fusing info on BAMs and VCF files so we know which VCF file goes with which BAM files and BAM object ids
# For info on downloaded data see notes-on-icgc.txt


bam.info <- data.table::fread("data-raw/collaboratory_WGS_BAMS_repository_1623246282.tsv")

stopifnot(unique(bam.info$format) == "BAM")



bam.info <- bam.info[  ,
                       c("Object ID", "File Name", "ICGC Donor", "Specimen ID", "Specimen Type",
                         "Size (bytes)")]

si <- split(bam.info, f = bam.info$`ICGC Donor`)

# Donors with multiple tumor samples

si[which(unlist(lapply(si, nrow)) > 4)]

# Need to drive by SP99999 ID from VCF?

organize.one.donor <- function(ssi) {

  by.spec <- split(ssi, ssi$`Specimen ID`)

  max.specs <- lapply(by.spec, function(one.spec) {
     max.idx <- which(one.spec$`Size (bytes)` == max(one.spec$`Size (bytes)`))
     if (length(max.idx) > 1) stop("max not unique")
     return(one.spec[max.idx, ])
  })

  # Find the normal specimens, not sure how to link to VCF file -- wait

  return(max.specs)
}

# debug(organize.one.donor)

donor.pairs <- lapply(si, organize.one.donor)

is.normal <- function(one.specimen) {
  return(grepl("Normal", one.specimen$`Specimen Type`))
}

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

merge.tumor.and.normal <- function(tumor.row, normal.row) {
  colnames(tumor.row) <- paste0("T_", colnames(tumor.row))
  colnames(normal.row) <- paste0("N_", colnames(normal.row))
  rr <- cbind(tumor.row, normal.row)
  return(rr)
}

# For debugging:
which(unlist(lapply(donor.pairs, length)) > 2)
# 302, 631, ...

no.match <- which(unlist(lapply(donor.pairs, length)) == 1)

donor.pairs <- donor.pairs[-no.match]

tumor.normal.pair <- tibble::as_tibble(do.call(rbind, lapply(donor.pairs, flatten.one.donor)))

samp.sheet <- data.table::fread("data-raw/pcawg_sample_sheet.v1.4.2016-09-14.tsv")

samp.sheet <- samp.sheet[samp.sheet$library_strategy == "WGS", ]

jj <- dplyr::inner_join(samp.sheet, tumor.normal.pair, by  = c("icgc_specimen_id" = "T_Specimen ID"))
# xjj <- dplyr::inner_join(samp.sheet, tumor.normal.pair, by  = c("icgc_specimen_id" = "N_Specimen ID"))

data.table::fwrite(jj, "collaboratory_bams_full.csv")

## Simplify the columns
jjj <- jj[ , c("aliquot_id", "icgc_donor_id", "icgc_specimen_id", "T_Object ID", "T_File Name", "N_Specimen ID", "N_Object ID", "N_File Name")]
colnames(jjj)[3] <- "T_Specimen ID"

data.table::fwrite(jjj, "collaboratory_bams.csv")
