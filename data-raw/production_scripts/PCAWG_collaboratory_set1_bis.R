#!/usr/local/bin/Rscript

library(DBSverify)

PCAWG_read_table_and_evaluate_DBS(
  in.table = "~/DBSverify/data-raw/production_scripts/collaboratory_bams_2021_07_16.csv",
  in.vcf.dir  = "~/mvv/PCAWG_all_Collaboratory_VCFs_and_BEDs/",
  minibam.dir = "~/mvv/collab-minibams-set1/",
  out.vcf.dir = "~/mvv/collab-minibams-set1-output-VCFs-bis2/")

# Just in case options(warn = 1) didn't work
warnings()
