#' Compare the contents of two directories of .sam files
#'
#' @param old.dir One directory.
#'
#' @param new.dir The other directory.
#'
#' For use as part of automated tests.
#'
#' @keywords internal

RegressSAMDirectory<- function(old.dir, new.dir) {
  new.sam.names <- dir(new.dir, pattern = ".sam")
  old.sam.names <- dir(old.dir, pattern = ".sam")
  testthat::expect_equal(new.sam.names, old.sam.names)
  for (ss in new.sam.names) {
    new <- ReadSamfile(file.path(new.dir, ss))
    old <- ReadSamfile(file.path(old.dir, ss))
    testthat::expect_equal(old, new)
  }
}
