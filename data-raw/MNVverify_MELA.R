## script to evaluate DBS calls in raw data

## workflow:
## Start with an ICAMS DBS vcf
## for each variant (row):
## subset the tumor and normal bam file for the region around the variant
## convert bam to sam
## remove reads that end half-way in the DBS (regardless of whether they are mutant or WT)
## for each remaining read, check whether they support the WT, DBS or half the DBS
## based on those counts, infer whether or not this event was likely to be a true DBS

########################################## setup ##########################################
setwd("../mvv")
## check whether there is a temp folder for intermediate files, if not, create

## specify which sample and bam-files to look at
smp <- "MELA-0161"
vcf.name <- "22.vcf"
Nbam <- "MELA-0161N_ch22.bam"
Tbam <- "MELA-0161T_ch22.bam"

##################################### helper functions #####################################

myPCAWGVAF<-function(vcf){
  tmp<-vcf[[8]]
  tmp<-stringr::str_split_fixed(tmp,"VAF=",2)[,2]
  tmp<-stringr::str_split_fixed(tmp,";",2)[,1]
  tmp<-as.numeric(as.character(tmp))
  return(cbind(vcf, VAF=tmp, read.depth = rep(0,length(tmp))))
}

vcf.list <-ICAMS::ReadAndSplitVCFs(vcf.name,
                 variant.caller = "unknown",
                 num.of.cores = 10,
                 names.of.VCFs = smp,
                 get.vaf.function = myPCAWGVAF,
                 max.vaf.diff = 1)

DBS.vcf <- vcf.list$DBS[[1]]

evaluated.vcf <- VerifyDBSVcf(DBS.vcf, Nbam = Nbam, Tbam = Tbam)

outfile <- paste0(smp,"_PCAWG_DBS_evaluated.vcf")
write.table(evaluated.vcf, file = outfile, sep="\t", quote=F, row.names = F)

old <- data.table::fread("MELA-0161_PCAWG_DBS_evaluated.vcf")
new <- data.table::fread(outfile)
all.equal(old, new)
cat("regression result", all.equal(old, new), "\n")

