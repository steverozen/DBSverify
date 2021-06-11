files = c(
  "~/collab.vcf/0009b464-b376-4fbc-8a56-da538269a02f.broad-mutect-v3.20160222.somatic.snv_mnv.vcf.gz",
  "~/collab.vcf/0009b464-b376-4fbc-8a56-da538269a02f.dkfz-snvCalling_1-0-132-1-hpc.1510221050.somatic.snv_mnv.vcf.gz",
  "~/collab.vcf/0009b464-b376-4fbc-8a56-da538269a02f.MUSE_1-0rc-vcf.20151216.somatic.snv_mnv.vcf.gz",
  "~/collab.vcf/0009b464-b376-4fbc-8a56-da538269a02f.svcp_1-0-4.20150204.somatic.snv_mnv.vcf.gz",
  "~/pcawg.vcf/final_consensus_12aug_passonly/snv_mnv//0009b464-b376-4fbc-8a56-da538269a02f.consensus.snv_mnv.vcf.gz"


)

library(dplyr)
library(tibble)

dbs.foo <- ICAMS::ReadAndSplitVCFs(files = files, variant.caller = "unknown", always.merge.SBS = TRUE )

pfoo <- dbs.foo$DBS[[5]]
which(pfoo$CHROM == 15 & pfoo$POS == 40527158)

pfoo <- tibble(dbs.foo$DBS[[5]])
which(pfoo$CHROM == 15 & pfoo$POS == 40527158)

dbs <- list(mutect=tibble(dbs.foo$DBS[[1]][ , 1:4]),
            dkfz  =tibble(dbs.foo$DBS[[2]][ , 1:4]),
            muse  =tibble(dbs.foo$DBS[[3]][ , 1:4]),
            sanger=tibble(dbs.foo$DBS[[4]][ , 1:4]),
            pcawg =tibble(dbs.foo$DBS[[5]]))

pfoo <- dbs$pcawg
which(pfoo$CHROM == 15 & pfoo$POS == 40527158)
which(dbs$pcawg$CHROM == 15 & dbs$pcawg$POS == 40527158)

colnames(dbs$mutect)[3:4] <- paste0(colnames(dbs$mutect)[3:4], "_mt")
colnames(dbs$dkfz)[3:4] <- paste0(colnames(dbs$dkfz)[3:4], "_dk")
colnames(dbs$muse)[3:4] <- paste0(colnames(dbs$muse)[3:4], "_ms")
colnames(dbs$sanger)[3:4] <- paste0(colnames(dbs$sanger)[3:4], "_sa")

all.dbs <- full_join(dbs$mutect,
                     full_join(dbs$dkfz,
                               full_join(dbs$muse,
                                         full_join(dbs$sanger, dbs$pcawg))))

one.row <- function(x) {
  return(sum(!is.na(x[3:10]))/2)
}

# debug(one.row)

dnc <- function(df, col) {
  ii <- which(colnames(df) == col)
  if (length(ii) > 0) {
    df <- df[ , -ii]
  }
  return(df)
}

all.dbs <- dnc(all.dbs, "VAF")
all.dbs <- dnc(all.dbs, "read.depth")
all.dbs <- dnc(all.dbs, "remark.for.DBS")
all.dbs <- dnc(all.dbs, "ID")
all.dbs <- dnc(all.dbs, "INFO")
all.dbs <- dnc(all.dbs, "FORMAT")
all.dbs <- dnc(all.dbs, "AOCS-117-9-AOCS-117-13")


all.dbs$num.support <- apply(all.dbs, FUN = one.row, MARGIN = 1)

dim(all.dbs)
table(all.dbs$num.support, !is.na(all.dbs$REF), dnn = c("num.caller.supporting", "called.by.pcawg"))


sbs.foo <- ICAMS::ReadAndSplitVCFs(files = files, variant.caller = "unknown", always.merge.SBS = FALSE )



sbs <- list(mutect=tibble(sbs.foo$SBS[[1]][ , c(1,2,4,5)]),
            dkfz  =tibble(sbs.foo$SBS[[2]][ , c(1,2,4,5)]),
            muse  =tibble(sbs.foo$SBS[[3]][ , c(1,2,4,5)]),
            sanger=tibble(sbs.foo$SBS[[4]][ , c(1,2,4,5)]),
            pcawg =tibble(sbs.foo$SBS[[5]][ , c(1,2,4:8)]))

colnames(sbs$mutect)[3:4] <- paste0(colnames(sbs$mutect)[3:4], "_SBS_mt")
colnames(sbs$dkfz)[3:4] <- paste0(colnames(sbs$dkfz)[3:4], "_SBS_dk")
colnames(sbs$muse)[3:4] <- paste0(colnames(sbs$muse)[3:4], "_SBS_ms")
colnames(sbs$sanger)[3:4] <- paste0(colnames(sbs$sanger)[3:4], "_SBS_sa")
colnames(sbs$pcawg)[3:ncol(sbs$pcawg)] <- paste0(colnames(sbs$pcawg)[3:ncol(sbs$pcawg)], "_SBS_pc")

all.sbs <- full_join(sbs$mutect,
                     full_join(sbs$dkfz,
                               full_join(sbs$muse,
                                         full_join(sbs$sanger, sbs$pcawg))))

all <- full_join(all.dbs, all.sbs)

alls <- (all %>% arrange(CHROM, POS))

dbsidx <- which(!is.na(alls$num.support))

dd <- c(dbsidx, (dbsidx + 1))
dd <- sort(dd)

alldd <- (alls[dd, ])

idx2 <- which(!is.na(alldd$num.support) & alldd$num.support > 1 & is.na(alldd$REF)) # NOT IN PCAWG BUT OTHER CALLERS HAVE
idd <- c(idx2, idx2 + 1)
idd <- sort(idd)

investigate <- alldd[idd, ]
View(investigate[53:54, ])  # why did't the PCAWG SBSs get merged? # They did get merged, see below
# 15 : 40527158
# 15 : 40527159

pfoo <- sbs.foo$SBS[[5]]
which(pfoo$CHROM == 15 & pfoo$POS == 40527158)
which(pfoo$CHROM == 15 & pfoo$POS == 40527159)


pfoo <- dbs.foo$SBS[[5]]
which(pfoo$CHROM == 15 & pfoo$POS == 40527158)
which(pfoo$CHROM == 15 & pfoo$POS == 40527159)

pfoo <- dbs.foo$DBS[[5]]
which(pfoo$CHROM == 15 & pfoo$POS == 40527158)
which(pfoo$CHROM == 15 & pfoo$POS == 40527159)


