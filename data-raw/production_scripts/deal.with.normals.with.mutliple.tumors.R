
# move bad normal minibams to a subdir
# (bad because they are minibams for multiple tumor samples, but don't have
# all the DBS sites for all the tumors)
f1 <- function(path1 = "~/mvv/collab-minibams-set1/") {
  bb <- data.table::fread("~/DBSverify/data-raw/production_scripts/multitumor.DBSverify.info.2021.08.03.csv")
  bb <- bb[ , c(-1, -3)]
  path2 <- file.path(path1, "bad.normal.minibams")
  for (ii in 1:nrow(bb)) {
    do <- bb$icgc.donor.id[ii]
    n.sp <- bb$N.specimen.id[ii]
    bam <- paste0(do, "_", n.sp, "_dbs_srt.bam")
    bami <- paste0(bam, ".bai")
    bam.path <- file.path(path1, bam)
    if (file.exists(bam.path)) {
      cat(bam.path, "\n")
      system2("mv", c(bam.path, path2))
      system2("mv", c(file.path(path1, bami), path2))
    }
  }
}


f1()
f1("~/mvv/collab-minibams-set2/")
f1("~/mvv/collab-minibams-set3/")


f2 <- function(path1 = "~/mvv/collab-minibams-set1/") {
  bb <- data.table::fread("~/DBSverify/data-raw/production_scripts/multitumor.DBSverify.info.2021.08.03.csv")
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

f2()
f2("~/mvv/collab-minibams-set2/")
f2("~/mvv/collab-minibams-set3/")
