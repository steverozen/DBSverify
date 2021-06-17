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


VCF_to_BED <- function(in.vcf, out.bed, padding = 10) {
  stopifnot(data.table::is.data.table(in.vcf) || tibble::is_tibble(in.vcf))
  bed.table <- in.vcf[ , 1:2, drop = FALSE]
  bed.table$start <- in.vcf[ , 2] - (padding + 1)
  bed.table$end   <- in.vcf[ , 2] + padding
  bed.table <- bed.table[ , -2]
  data.table::fwrite(data.table::as.data.table(bed.table), out.bed, sep = "\t", col.names = FALSE)
  return(invisible(bed.table))
}


## Stuff to give to Willie on 14 Jun
one.sup.dbs  <- VCF_to_BED(all.dbs[all.dbs$num.support == 1, ], "one.support.dbs.bed")
two.sup.dbs  <- VCF_to_BED(all.dbs[all.dbs$num.support == 2, ], "two.support.dbs.bed")
high.sup.dbs <- VCF_to_BED(all.dbs[all.dbs$num.support > 2, ],  "high.support.dbs.bed")

# test.case <- which(samp.sheet$aliquot_id == "0009b464-b376-4fbc-8a56-da538269a02f")

# samp.sheet[test.case, ]$icgc_donor_id
#"DO46416"

# from other .R file

donor.pairs[["DO46416"]]

# $SP101724
# Object ID                            File Name ICGC Donor Specimen ID            Specimen Type Size (bytes)
# 1: 2c9d90f1-5c46-5f3f-8a7f-771b36414cf2 bf5b874240ec430239013f4292d4151c.bam    DO46416    SP101724 Recurrent tumour - other 117118114996

# $SP101728
# Object ID                            File Name ICGC Donor Specimen ID          Specimen Type Size (bytes)
# 1: 9ec1f43b-379c-58cd-85fe-3d09439058f1 7441677a3aed0864a23b17b84c0a28c9.bam    DO46416    SP101728 Normal - blood derived 134788648263


merge.callers <- function(row) {
  row <- unlist(row)
  ref <- row[c("REF_mt", "REF_dk", "REF_ms", "REF_sa")]
  alt <- row[c("ALT_mt", "ALT_dk", "ALT_ms", "ALT_sa")]
  ref.ok <- which(!is.na(ref))
  alt.ok <- which(!is.na(alt))
  ref2 <- unique(ref[ref.ok])
  alt2 <- unique(alt[alt.ok])
  stopifnot(length(ref2) == 1)
  stopifnot(length(alt2) == 1)
  rr <- c(row[c("CHROM", "POS")], REF = ref2, ALT = alt2, num.support = length(ref.ok), pcawg.call = !is.na(row["REF"]))
  return(rr)
}


## Test 3 kinds of unioned DBSs -- only one caller supported, two caller supported, and 3 or 4 callers supported
setwd("~/mvv/test.minibams/")

tmp.vcf <- all.dbs[all.dbs$num.support == 1, ]
# tmp.vcf <- tmp.vcf[1:10, ]
tmp2.vcf <- t(apply(tmp.vcf, MARGIN = 1, FUN = merge.callers))
colnames(tmp2.vcf) <- c("#CHROM", "POS", "REF", "ALT", "num.callers", "pcawg.called")
data.table::fwrite(tmp2.vcf, "tmp.vcf", sep = "\t")
new.one.sup <-
  DBSverify::Read_DBS_VCF_and_BAMs_to_verify_DBSs(
  input.vcf = "tmp.vcf",
  Nbam.name = "SP101728_oneSupport.bam",
  Tbam.name = "SP101724_oneSupport.bam",
  N.slice.dir = "one.support.N.slice.dir",
  T.slice.dir = "one.support.T.slice.dir",
  unlink.slice.dir = FALSE,
  verbose = 1,
  outfile = "new.one.support.eval.vcf"
)

# Test with something like
# fisher.test(matrix(c(40,0,35,3), ncol = 2), alternative = "g")

tmp.vcf <- all.dbs[all.dbs$num.support == 2, ]
tmp2.vcf <- t(apply(tmp.vcf, MARGIN = 1, FUN = merge.callers))
colnames(tmp2.vcf) <- c("#CHROM", "POS", "REF", "ALT", "num.callers", "called.by.pcawg")
data.table::fwrite(tmp2.vcf, "two.tmp.vcf", sep = "\t")
two.sup <-
  DBSverify::Read_DBS_VCF_and_BAMs_to_verify_DBSs(
    input.vcf = "two.tmp.vcf",
    Nbam.name = "SP101728_twoSupport.bam",
    Tbam.name = "SP101724_twoSupport.bam",
    N.slice.dir = "two.support.N.slice.dir",
    T.slice.dir = "two.support.T.slice.dir",
    unlink.slice.dir = FALSE,
    verbose = 1,
    outfile = "old.two.support.eval.vcf"
  )

# Investigate? 8:126314418
# " 4 67785957

tmp.vcf <- all.dbs[all.dbs$num.support > 2, ]
tmp2.vcf <- t(apply(tmp.vcf, MARGIN = 1, FUN = merge.callers))
colnames(tmp2.vcf) <- c("#CHROM", "POS", "REF", "ALT", "num.callers", "called.by.pcawg")
data.table::fwrite(tmp2.vcf, "high.tmp.vcf", sep = "\t")
high.sup <-
  DBSverify::Read_DBS_VCF_and_BAMs_to_verify_DBSs(
    input.vcf = "high.tmp.vcf",
    Nbam.name = "SP101728_highSupport.bam",
    Tbam.name = "SP101724_highSupport.bam",
    N.slice.dir = "high.support.N.slice.dir",
    T.slice.dir = "high.support.T.slice.dir",
    unlink.slice.dir = FALSE,
    verbose = 1,
    outfile = "old.high.support.eval.vcf"
  )

setwd("~/DBSverify/")


## From here down is an investigation of SBSs calls and their relationship to DBS calls.
## Very few odd cases, so we just went with DBSs calls (as synthesized from individual caller SBS calls; i.e.
## for each caller we generated DBS calls, and then proceeded from there.)

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

sup.2.and.pcawg <- which(is.na(alls$REF) & alls$num.support == 2)
sup.2.and.pcawg<- sort(c(sup.2.and.pcawg + 1, sup.2.and.pcawg))

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
