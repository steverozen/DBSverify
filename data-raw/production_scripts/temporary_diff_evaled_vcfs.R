setwd("~/mvv/collab-minibams-set1-output-VCFs-bis2/")
vcfs <- dir(pattern = "\\.vcf")
for (vv in vcfs) {
  xx <- gsub("SP\\d+_", "", vv, perl = TRUE)
   system2("diff",
           c(file.path("~/mvv/collab-minibams-set1-output-VCFs-bis/", xx),
             vv))
}


set1 <- dir("~/mvv/collab-minibams-set1/", pattern = "\\.bam")
set2<- dir("~/mvv/collab-minibams-set2/", pattern = "\\.bam")
setdiff(set1, set2)
only.set2 <- setdiff(set2, set1)
dups <- intersect(set1, set2)

for (ss in dups) {
  system2("mv",
          c(file.path("~/mvv/collab-minibams-set2/", ss),
              "~/mvv/collab-minibams-set2/dups.w.set1/"))
}

# move bad normal minibams to a subdir
# (bad because they are minibams for multiple tumor samples, but don't have
# all the DBS sites for all the tumors)
f1 <- function() {
  bb <- data.table::fread("~/DBSverify/data-raw/production_scripts/multitumor.DBSverify.info.2021.08.03.csv")
  bb <- bb[ , c(-1, -3)]
  path1 <- "~/mvv/collab-minibams-set1/"
  path2 <- file.path(path1, "bad.normal.minibams")
  for (ii in 1:nrow(bb)) {
    do <- bb$icgc.donor.id[ii]
    n.sp <- bb$N.specimen.id[ii]
    bam <- paste0(do, "_", n.sp, "_dbs_srt.bam")
    bami <- paste0(bam, ".bai")
    bam.path <- file.path(path1, bam)
    if (file.exists(bam.path)) {
      cat(bam.path)
      system2("mv", c(bam.path, path2))
      system2("mv", c(file.path(path1, bami), path2))
    }
  }
}

f2 <- function() {
  bb <- data.table::fread("~/DBSverify/data-raw/production_scripts/multitumor.DBSverify.info.2021.08.03.csv")
  path1 <- path1 <- "~/mvv/collab-minibams-set2/"
  path2 <- file.path(path1, "redo.tumor.minibams")
  for (ii in 1:nrow(bb)) {
    do <- bb$icgc.donor.id[ii]
    t.sp <- bb$T.specimen.id[ii]
    bam <- paste0(do, "_", t.sp, "_dbs_srt.bam")
    bami <- paste0(bam, ".bai")
    bam.path <- file.path(path1, bam)
    if (file.exists(bam.path)) {
      cat(bam.path, "\n")
      system2("mv", c(bam.path, path2))
      system2("mv", c(file.path(path1, bami), path2))
    }

  }
}
