# Keep track of new versus already analyzed miniBAMs

sort.minibams <- function() {
  set1 <- dir("~/mvv/collab-minibams-set1/", pattern = "\\.bam")
  # Includes .bam and .bami files
  set2 <- dir("~/mvv/collab-minibams-set2/", pattern = "\\.bam")
  new.set.path <- "~/mvv/collab-minibams-set3/"
  set3 <- dir(new.set.path, pattern = "\\.bam")
  prev.set <- c(set1, set2)
  stopifnot(!any(duplicated(prev.set)))
  only.new <- setdiff(set3, prev.set)
  dups <- intersect(prev.set, set3)

  for (ss in dups) {
    system2("mv",
            c(file.path(new.set.path, ss),
              file.path(new.set.path, "dups.w.prev")))
  }
}
debug(sort.minibams)
