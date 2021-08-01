library(ICAMS)
options(warn = 1)
vcf.dir <- "/home/gmssgr/mvv/collab-minibams-set1-output-VCFs/"
vcfs <- dir(vcf.dir, full.names = TRUE)
vcf.aliquot.id1 <- gsub(vcf.dir, "", vcfs)
vcf.aliquot.id  <- gsub("/DO\\d+_", "", vcf.aliquot.id1, perl = TRUE)
aliquot.id <- gsub("_PCAWG_evaluated.vcf", "", vcf.aliquot.id)
tt <- data.table::fread(
  "~/DBSverify/data-raw/production_scripts/collaboratory_bams_2021_07_16.csv")

jj <- dplyr::left_join(
  data.table::data.table(aliquot_id = aliquot.id), tt)


PCAWG_filter_DBS <- function(VCF.file) {
  vcf <- data.table::fread(VCF.file)
  colnames(vcf)[1] <- "CHROM"
  vcf$CHROM <- as.character(vcf$CHROM)
  orig.n.dbs <- nrow(vcf)
  pcawg.vcf <- vcf[grep("pcawg", vcf$which.callers), ]
  orig.pcawg.dbs <- nrow(pcawg.vcf)
  final.vcf <- pcawg.vcf[pcawg.vcf$DBSconclusion == "True DBS", ]
  return(list(final.vcf      = final.vcf,
              pcawg.vcf      = pcawg.vcf,
              orig.n.dbs     = orig.n.dbs,
              orig.pcawg.dbs = orig.pcawg.dbs,
              good.n.dbs     = nrow(vcf)))
}
# debug(PCAWG_filter_DBS)
xx <- lapply(vcfs, PCAWG_filter_DBS)

final.vcf.list <- lapply(xx, function(xx) xx$final.vcf)
names(final.vcf.list) <- jj$`T_Specimen ID`

dbs.cats <- ICAMS::VCFsToDBSCatalogs(
  final.vcf.list,
  ref.genome = "hg19",
)


ICAMS::PlotCatalogToPdf(dbs.cats$catDBS78, "collab-set1-DBS.pdf")

if (FALSE) {
  # ask nanhai to debug
all.output <- ICAMS::VCFsToCatalogs(
  files = vcfs,
  ref.genome = "hg19",
  num.of.cores = 1,
  region = "genome"
)
}

