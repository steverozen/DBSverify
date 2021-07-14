if (FALSE) {
  HMF_prep_BEDs_and_DBS_VCFs(input.vcf.dir = "~/mvv/HMF-VCFs/")
}

HMF_prep_BEDs_and_DBS_VCFs <- function(input.vcf.dir,
                                       out.bed.dir = input.vcf.dir,
                                       out.vcf.dir = input.vcf.dir,
                                       verbose = TRUE) {

  if (verbose) message("HMV VCFs from ", input.vcf.dir)

  input.vcfs <- dir(input.vcf.dir,
                    full.names = FALSE,
                    # pattern = "CPCT02030256T.purple.somatic.vcf.gz")
                    pattern = ".vcf.gz")

  # ignore possible index files

  input.vcfs <-
    grep("(\\.idx$)|(\\.tbi$)", input.vcfs, value = TRUE, invert = TRUE)

  for (in.vcf in input.vcfs) {
    HMF_prep_1_BED_and_DBS_VCF(in.vcf,
                               out.bed.dir = out.bed.dir,
                               out.vcf.dir = out.vcf.dir,
                               verbose     = verbose,
                               input.vcf.dir = input.vcf.dir)

  }
}


HMF_prep_1_BED_and_DBS_VCF <- function(in.vcf,
                                       out.bed.dir,
                                       out.vcf.dir,
                                       verbose,
                                       input.vcf.dir) {

  if (verbose) message("Processing ", in.vcf)

  vcf <- ICAMS:::ReadVCFs(file             = file.path(input.vcf.dir, in.vcf),
                         variant.caller   = "unknown",
                         filter.status    = "PASS")

  vcf <- vcf[[1]]
  vcf <- Remove_non_canonical_chromosomes(vcf)

  # Discard indels
  vcf <- vcf[-(which(nchar(vcf$REF) != nchar(vcf$ALT))), ]
  if (nrow(vcf) == 0) {
    if (verbose) message("No rows in in.vcf")
    return(NULL)
  }
  vcf <- remove_named_column(vcf, c("FILTER", "FORMAT", "QUAL"))

  dbs <- vcf[which(nchar(vcf$REF) == 2), ]
  if (nrow(dbs) > 0) {
    dbs <- dbs[ , c("CHROM", "POS", "REF", "ALT")]
    dbs$source <- "original"
    if (verbose) message(nrow(dbs), " Original DBSs")
  } else {
    dbs <- NULL
  }
  sbs <- vcf[which(nchar(vcf$REF) == 1), ]
  if (nrow(sbs) > 0 ) {
    vvv <- ICAMS:::SplitSBSVCF(sbs, always.merge.SBS = TRUE)
    dbs.from.sbs <- vvv$DBS.vcf
    if (nrow(dbs.from.sbs > 0)) {
      dbs.from.sbs <- dbs.from.sbs[ , c("CHROM", "POS", "REF", "ALT")]
      if (verbose) message(nrow(dbs.from.sbs), " DBSs from SBS")
      dbs.from.sbs$source <- "adjacent_SBS"
    } else {
      dbs.from.sbs <- NULL
    }
  } else {
    dbs.from.sbs <- NULL
  }
  if (!is.null(dbs) & !is.null(dbs.from.sbs)) {
    out.dbs <- rbind(dbs, dbs.from.sbs)
  } else if (!is.null(dbs)) {
    out.dbs <- dbs
  } else {
    out.dbs <- dbs.from.sbs
  }

  file.root <- gsub(".purple.somatic.vcf.gz", "", in.vcf, fixed = TRUE)
  if (verbose) message("Writing bed file")
  bed <- VCF_to_BED(
    out.dbs,
    out.bed = file.path(out.bed.dir,
                        paste0(file.root, "_HMF_DBS.bed")))

  if (verbose) {
    message("Writing vcf file: ", nrow(out.dbs), " rows")
  }

  data.table::fwrite(out.dbs,
                     file.path(out.vcf.dir,
                               paste0(file.root, "_HMF_DBS.vcf")),
                     sep = "\t")


  return(invisible(out.dbs))
}
