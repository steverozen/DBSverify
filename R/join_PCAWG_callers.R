join_PCAWG_callers <- function(aliquot.id, indiv.vcf.dir, pcawg.vcf.dir) {


  input.vcfs <- dir(indiv.vcf.dir,
                    pattern = "0009b464-b376-4fbc-8a56-da538269a02f",
                    full.names = TRUE)
  if (length(input.vcfs) != 8) {
    stop("Not enough input VCFs for aliquot id", aliquot.id)
  }

  br <- grep("(\\.idx$)|(\\.tbi$)",
             grep("\\.broad", input.vcfs, value = TRUE), value = TRUE, invert = TRUE)
  dk <- grep("(\\.idx$)|(\\.tbi$)",
             grep("\\.dkfz", input.vcfs, value = TRUE), value = TRUE, invert = TRUE)
  mu <- grep("(\\.idx$)|(\\.tbi$)",
             grep("\\.MUSE", input.vcfs, value = TRUE), value = TRUE, invert = TRUE)
  sv <- grep("(\\.idx$)|(\\.tbi$)",
             grep("\\.svcp", input.vcfs, value = TRUE), value = TRUE, invert = TRUE)

  pc <-     file.path(pcawg.vcf.dir,
                      paste0(aliquot.id, ".consensus.snv_mnv.vcf.gz"))

  dbs.foo <- ICAMS::ReadAndSplitVCFs(files = c(br, dk, mu, sv, pc), variant.caller = "unknown", always.merge.SBS = TRUE )

  dbs <- list(mutect = tibble::tibble(dbs.foo$DBS[[1]][ , 1:4]),
              dkfz   = tibble::tibble(dbs.foo$DBS[[2]][ , 1:4]),
              muse   = tibble::tibble(dbs.foo$DBS[[3]][ , 1:4]),
              sanger = tibble::tibble(dbs.foo$DBS[[4]][ , 1:4]),
              pcawg  = tibble::tibble(dbs.foo$DBS[[5]]))

  # Mark the column names so we know which column came from which caller
  colnames(dbs$mutect)[3:4] <- paste0(colnames(dbs$mutect)[3:4], "_mt")
  colnames(dbs$dkfz)[3:4] <- paste0(colnames(dbs$dkfz)[3:4], "_dk")
  colnames(dbs$muse)[3:4] <- paste0(colnames(dbs$muse)[3:4], "_ms")
  colnames(dbs$sanger)[3:4] <- paste0(colnames(dbs$sanger)[3:4], "_sa")

  all.dbs <- dplyr::full_join(
    dbs$mutect,
    dplyr::full_join(
      dbs$dkfz,
      dplyr::full_join(
        dbs$muse,
        dplyr::full_join(dbs$sanger, dbs$pcawg))))

  one.row <- function(x) {
    return(sum(!is.na(x[3:10]))/2)
  }

  # Function to delete columns by name
  dnc <- function(df, col) {
    ii <- which(colnames(df) == col)
    if (length(ii) > 0) {
      df <- df[ , -ii]
    }
    return(df)
  }

  all.dbs <- dnc(all.dbs, "VAF")
  all.dbs <- dnc(all.dbs, "read.depth")
  all.dbs <- dnc(all.dbs, "remark.for.DBS")
  all.dbs <- dnc(all.dbs, "ID")
  all.dbs <- dnc(all.dbs, "INFO")
  all.dbs <- dnc(all.dbs, "FORMAT")
  all.dbs <- dnc(all.dbs, "AOCS-117-9-AOCS-117-13")

  all.dbs$num.support <- apply(all.dbs, FUN = one.row, MARGIN = 1)




  return(all.dbs)
}
