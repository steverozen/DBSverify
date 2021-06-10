# Join info on BAM file object ids with SP IDs

samp.sheet <- data.table::fread("data-raw/pcawg_sample_sheet.v1.4.2016-09-14.tsv")

bam.info <- data.table::fread("data-raw/collaboratory_WGS_BAMS_repository_1623246282.tsv")

stopifnot(unique(bam.info$format) == "BAM")



bam.info <- bam.info[  ,
                       c("Object ID", "ICGC Donor", "Specimen ID", "Specimen Type",
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

debug(organize.one.donor)

foo <- lapply(si, organize.one.donor)
