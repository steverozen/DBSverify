setwd("~/mvv/collab-minibams-set1-output-VCFs-bis2/")
vcfs <- dir(pattern = "\\.vcf")
for (vv in vcfs) {
  xx <- gsub("SP\\d+_", "", vv, perl = TRUE)
   system2("diff",
           c(file.path("~/mvv/collab-minibams-set1-output-VCFs-bis/", xx),
             vv))
}
