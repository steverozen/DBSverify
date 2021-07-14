PCAWG_join_callers <- function(aliquot.id, indiv.vcf.dir, pcawg.vcf.dir, verbose = TRUE) {

  if (verbose) message("aliquot.id = ", aliquot.id)

  input.vcfs <- dir(indiv.vcf.dir,
                    pattern = aliquot.id,
                    full.names = TRUE)

  # ignore possible index files

  input.vcfs <-
    grep("(\\.idx$)|(\\.tbi$)", input.vcfs, value = TRUE, invert = TRUE)

  if (length(input.vcfs) != 4) {
    message("Not enough input VCFs for aliquot id", aliquot.id)
    # Data frame with 0 rows
    return(data.frame(1)[-1, , drop = FALSE])
  }

  br <- grep("\\.broad", input.vcfs, value = TRUE)
  dk <- grep("\\.dkfz", input.vcfs, value = TRUE)
  mu <- grep("\\.MUSE", input.vcfs, value = TRUE)
  sv <- grep("\\.svcp", input.vcfs, value = TRUE)

  pc.path <- dir(pcawg.vcf.dir,
                 pattern = aliquot.id,
                 full.names = TRUE)

  pc <- grep("(\\.idx$)|(\\.tbi$)", pc.path, value = TRUE, invert = TRUE)

  if (length(pc) != 1) {
    stop("Cannot find PCAWG consensus VCF file")
  }

  pcawg.all <- ICAMS::ReadAndSplitVCFs(files            = pc,
                                       variant.caller   = "unknown",
                                       filter.status    = NULL,
                                       always.merge.SBS = TRUE)
  pcawg.vcf  <- Remove_non_canonical_chromosomes(pcawg.all$DBS[[1]])
  rm(pc, pcawg.all)
  if (verbose) message("nrow(pcawg.vcf) = ", nrow(pcawg.vcf))

  other.all <- ICAMS::ReadAndSplitVCFs(files            = c(br, dk, mu, sv),
                                       variant.caller   = "unknown",
                                       filter.status    = "PASS",
                                       always.merge.SBS = TRUE )

  mutect.vcf <- Remove_non_canonical_chromosomes(other.all$DBS[[1]][ , 1:4])
  if (verbose) message("nrow(mutect.vcf) = ", nrow(mutect.vcf))
  dkfz.vcf   <- Remove_non_canonical_chromosomes(other.all$DBS[[2]][ , 1:4])
  if (verbose) message("nrow(dkvz.vcf) = ", nrow(dkfz.vcf))
  muse.vcf   <- Remove_non_canonical_chromosomes(other.all$DBS[[3]][ , 1:4])
  sanger.vcf <- Remove_non_canonical_chromosomes(other.all$DBS[[4]][ , 1:4])
  rm(other.all, br, dk, mu, sv)

  dbs <- list(mutect = tibble::tibble(mutect.vcf),
              dkfz   = tibble::tibble(dkfz.vcf),
              muse   = tibble::tibble(muse.vcf),
              sanger = tibble::tibble(sanger.vcf),
              pcawg  = tibble::tibble(pcawg.vcf))

  # Mark the column names so we know which column came from which caller
  colnames(dbs$mutect)[3:4] <- paste0(colnames(dbs$mutect)[3:4], "_mt")
  colnames(dbs$dkfz)[3:4]   <- paste0(colnames(dbs$dkfz)[3:4], "_dk")
  colnames(dbs$muse)[3:4]   <- paste0(colnames(dbs$muse)[3:4], "_ms")
  colnames(dbs$sanger)[3:4] <- paste0(colnames(dbs$sanger)[3:4], "_sa")
  colnames(dbs$pcawg)[3:4] <- paste0(colnames(dbs$pcawg[3:4]), "_pcawg")

  all.dbs <- dplyr::full_join(
    dbs$mutect,
    dplyr::full_join(
      dbs$dkfz,
      dplyr::full_join(
        dbs$muse,
        dplyr::full_join(dbs$sanger, dbs$pcawg))))

  #one.row <- function(x) {
  #  return(sum(!is.na(x[3:10]))/2)
  #}

  all.dbs <- remove_named_column(
    all.dbs,
    c("VAF", "read.depth", "remark.for.DBS", "ID", "INFO", "FORMAT",
      "AOCS-117-9-AOCS-117-13"))

  # all.dbs$num.support <- apply(all.dbs, FUN = one.row, MARGIN = 1)

  return(all.dbs)
}
