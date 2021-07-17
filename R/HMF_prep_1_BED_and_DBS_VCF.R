

HMF_prep_1_BED_and_DBS_VCF <- function(in.vcf,
                                       out.bed.dir,
                                       out.vcf.dir,
                                       verbose,
                                       input.vcf.dir) {

  if (verbose) message("Processing ", in.vcf)

  num.orignal.dbs  <- 0
  num.adjacent.sbs <- 0

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
    num.orignal.dbs <- nrow(dbs)
    if (verbose) message(num.orignal.dbs, " Original DBSs")

  } else {
    dbs <- NULL
  }
  sbs <- vcf[which(nchar(vcf$REF) == 1), ]
  if (nrow(sbs) > 0 ) {
    vvv <- ICAMS:::SplitSBSVCF(sbs, always.merge.SBS = TRUE)
    dbs.from.sbs <- vvv$DBS.vcf
    if (nrow(dbs.from.sbs > 0)) {
      dbs.from.sbs <- dbs.from.sbs[ , c("CHROM", "POS", "REF", "ALT")]
      num.adjacent.sbs <- nrow(dbs.from.sbs)
      if (verbose) message(num.adjacent.sbs, " DBSs from SBS")
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


  return(data.frame(in.vcf = in.vcf,
                    num.orignal.dbs = num.orignal.dbs,
                    num.adjacent.sbs = num.adjacent.sbs))
}
