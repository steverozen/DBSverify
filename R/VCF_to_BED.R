VCF_to_BED <- function(in.vcf, out.bed, padding = 10) {
  stopifnot(data.table::is.data.table(in.vcf) || tibble::is_tibble(in.vcf))
  in.vcf <- Remove_non_canonical_chromosomes(in.vcf)
  bed.table <- in.vcf[ , 1:2, drop = FALSE]
  bed.table$start <- in.vcf[ , 2] - (padding + 1)
  bed.table$end   <- in.vcf[ , 2] + padding
  bed.table <- bed.table[ , -2]
  data.table::fwrite(data.table::as.data.table(bed.table),
                     out.bed,
                     sep       = "\t",
                     col.names = FALSE,
                     scipen    = 10000000) # Don't use scientific notation
  return(invisible(bed.table))
}
