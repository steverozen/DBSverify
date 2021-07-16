beds <- dir("~/mvv/PCAWG_all_Collaboratory_VCFs_and_BEDs/", pattern = "\\.bed$")
tt <- data.table::fread("~/DBSverify/data-raw/collaboratory_bams_2021_07_13.csv")

bbb <- gsub("_merged_PCAWG_DBS.bed", "", beds)
bb <- tibble::tibble(aliquot_id = bbb)
xx <- dplyr::left_join(bb, tt)
