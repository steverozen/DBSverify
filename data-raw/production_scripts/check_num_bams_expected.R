# Check number of bams expected frmo output vcf

check.num.bams <- function(my.dir) {
    vcf <- dir(my.dir, pattern = ".vcf")
    xx <- strsplit(vcf, "_")
    yy <- t(matrix(unlist(xx), ncol = length(vcf)))
    yyy <- split(yy, yy[ , 1])

    num.normal.bams <- length(yyy)
    num.tumor.bams <- nrow(yy)
    message("Num tumor bams = ", num.tumor.bams,
            "\nNum normal bams = ", num.normal.bams,
            "\nTotal = ", num.normal.bams + num.tumor.bams)

  }

check.num.bams("~/mvv/collab-minibams-set1-output-VCFs/")
