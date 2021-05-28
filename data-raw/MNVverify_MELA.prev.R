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
i<-list.dirs(".",recursive = F)
if(!"./temp" %in% i){
  system("mkdir temp")
  print("NB: the newly created temp directory doesn't have the correct settings yet")
}

## specify which sample and bam-files to look at
smp <- "MELA-0161"
vcf.name <- "22.vcf"
Nbam <- "MELA-0161N_ch22.bam"
Tbam <- "MELA-0161T_ch22.bam"

##################################### helper functions #####################################

addID<-function(df,chr="CHROM",pos="POS",ref="REF",alt="ALT"){
  if(grepl("data.table",attributes(df)$class)){
    attributes(df)$class<-"data.frame"
  }
  df$ID<-paste0(df[,colnames(df) == chr],":",
                df[,colnames(df) == pos],":",
                df[,colnames(df) == ref],":",
                df[,colnames(df) == alt])
  return(df)
}

readSam<-function(file){
  df <- read.csv(file,sep="^",as.is=T,fill=T,header=F)
  df<-df[-grep("^@",df[,1]),]
  df<-as.data.frame(stringr::str_split_fixed(df,"\t",16))
  colnames(df)[1:11]<-c("QNAME","FLAG","CHROM","POS","MAPQ","CIGAR","Mate_CHROM","Mate_POS","InsertSize","SEQ","QUAL")

  ## remove weird rows if required
  df<-df[!is.na(suppressWarnings(as.numeric(df$FLAG))),]

  ## remove reads with flag > 256
  ## these are reads that
  ## - fail vendor QC
  ## - have multiple alignments
  ## - are marked as duplicates
  ## - or their pair aligns to a different chromosome
  df<-df[as.numeric(df$FLAG) < 256,]
  df<-df[df$Mate_CHROM == "=",]
  return(df)
}

# one variant
# returns a named character vector, the names are the read names, the values are .... plus left/right half supporting read
# Read supports only 2nd position"
# 1] "Read does not overlap complete DBS" "WT read"                            "Mut read"
#
readEval<-function(sam, v){
  tmp<-NULL
  for(r in 1:nrow(sam)){ # every row is a read
    ## check whether the entire DBS is covered by the read
    readRange<-c(as.numeric(sam$POS[r]),as.numeric(sam$POS[r])+nchar(sam$SEQ[r])-1)
    if(!(vcf[v,"POS"] >= min(readRange) & vcf[v,"POS"]+1 <= max(readRange))){
      tmp<-c(tmp,"Read does not overlap complete DBS")
      names(tmp)[length(tmp)]<-sam$QNAME[r]
    } else {

      ## only for script development checks
      ## sanity check; lookup BSgenome to confirm sequence matches (this should catch reads that should be inverted)
      #    Ranges <- GenomicRanges::GRanges(sam[r,"CHROM"],
      #                                     IRanges::IRanges(start = as.numeric(sam[r,"POS"]),
      #                                                      end = as.numeric(sam[r,"POS"]) + 20))
      #    checkSeq<-BSgenome::getSeq(ICAMS:::NormalizeGenomeArg("GRCh37"), Ranges, as.character = TRUE)
      #    substr(sam[r,"SEQ"],1,21) == checkSeq

      readSeq<-substr(sam[r,"SEQ"],vcf[v,"POS"]-as.numeric(sam[r,"POS"])+1,vcf[v,"POS"]-as.numeric(sam[r,"POS"])+1+1)

      ## check whether the sequence of the read matches the DBS
      pos1eval<-substr(readSeq,1,1) == substr(vcf[v,"REF"],1,1)
      pos2eval<-substr(readSeq,2,2) == substr(vcf[v,"REF"],2,2)
      if(pos1eval == T & pos2eval == T){
        tmp<-c(tmp,"WT read")
        names(tmp)[length(tmp)]<-sam$QNAME[r]
      } else {
        pos1eval<-substr(readSeq,1,1) == substr(vcf[v,"ALT"],1,1)
        pos2eval<-substr(readSeq,2,2) == substr(vcf[v,"ALT"],2,2)
        if(pos1eval == T & pos2eval == T){
          tmp<-c(tmp,"Mut read")
          names(tmp)[length(tmp)]<-sam$QNAME[r]
        } else {
          if(pos1eval == T & pos2eval == F){
            tmp<-c(tmp,"Read supports only 1st position")
            names(tmp)[length(tmp)]<-sam$QNAME[r]
          } else {
            if(pos1eval == F & pos2eval == T){
              tmp<-c(tmp,"Read supports only 2nd position")
              names(tmp)[length(tmp)]<-sam$QNAME[r]
            }
          }
        }
      }
    }
  }
  return(tmp)
}



myPCAWGVAF<-function(vcf){
  tmp<-vcf[[8]]
  tmp<-stringr::str_split_fixed(tmp,"VAF=",2)[,2]
  tmp<-stringr::str_split_fixed(tmp,";",2)[,1]
  tmp<-as.numeric(as.character(tmp))
  return(cbind(vcf, VAF=tmp, read.depth = rep(0,length(tmp))))
}

DBSconclusion<-function(vcf,germlineCutOff=0.2){
  vcf$DBSconclusion<-NA

  ## check whether one of the positions is actually a germline SNP
  tmp<-as.data.frame(stringr::str_split_fixed(vcf$NreadSupport,":",4))
  for(i in 1:4){tmp[,i]<-as.numeric(tmp[,i])}
  i<-sweep(tmp,1,rowSums(tmp),FUN="/")
  vcf$DBSconclusion[i[,2] >= germlineCutOff & i[,3] >= germlineCutOff]<-"DBS overlaps 2 germline SNPs"
  vcf$DBSconclusion[i[,2] >= germlineCutOff & i[,3] <= germlineCutOff]<-"1st position overlaps germline SNP"
  vcf$DBSconclusion[i[,2] <= germlineCutOff & i[,3] >= germlineCutOff]<-"2nd position overlaps germline SNP"
  vcf$DBSconclusion[i[,4] >= germlineCutOff ]<-"Germline DBS"
  ## for the evaluation of somatic events, don't look at rows that have germline involvement
  i<-is.na(vcf$DBSconclusion)

  ## Conclusions regarding somatic events
  tmp<-as.data.frame(stringr::str_split_fixed(vcf$TreadSupport,":",4))
  for(i in 1:4){tmp[,i]<-as.numeric(tmp[,i])}
  vcf$DBSconclusion[i & rowSums(tmp[,2:3]) <= 1 & tmp[,4] > 0 ]<-"True DBS"
  vcf$DBSconclusion[i & rowSums(tmp[,2:4]) == 0 ]<-"Neither position supported"
  vcf$DBSconclusion[i & rowSums(tmp[,2:3]) >1 ]<-"Adjacent SBSs"

  return(vcf)
}

######################################### load vcf #########################################

## load vcf
vcf.list <-ICAMS::ReadAndSplitVCFs(vcf.name,
                 variant.caller = "unknown",
                 num.of.cores = 10,
                 names.of.VCFs = smp,
                 get.vaf.function = myPCAWGVAF,
                 max.vaf.diff = 1)

write.table(vcf.list$DBS[[1]],file=paste0(smp,"_PCAWG_DBS.vcf"),sep="\t",row.names = F,quote=F)


## read DBS ICAMS vcf
vcf<-read.csv(paste0(smp,"_PCAWG_DBS.vcf"),sep="\t",as.is=T)

# vcf <- vcf.list$DBS[[1]]

#########
# vcf <- vcf[vcf$CHROM == "22", ]
#########

# debug(addID)
vcf<-addID(vcf)
## remove duplicates if present
vcf<-vcf[!duplicated(vcf$ID),]
rownames(vcf)<-vcf$ID

## prep columns for output
vcf$Format<-paste0("WtReads:pos1reads:pos2reads:MutReads")
vcf$NreadSupport<-NA
vcf$TreadSupport<-NA

################################### process each variant ###################################
tt <- function() {
for(v in vcf$ID){

  ## first create bam slices
  ## add a slight window around the variant, variants that don't overlap the (entire) DBS will be ignored anyway
  bamCoord <- paste0(vcf[v,"CHROM"],":",vcf[v,"POS"]-10,"-",vcf[v,"POS"]+10)
  system(paste0("samtools view -h ",Nbam," ",bamCoord," > temp/",smp,"_",bamCoord,"_N.sam"),wait=T)
  system(paste0("samtools view -h ",Tbam," ",bamCoord," > temp/",smp,"_",bamCoord,"_T.sam"),wait=T)

  ## process the sam files one by one
  for(s in c("T","N")){
    sam<-readSam(paste0("temp/",smp,"_",bamCoord,"_",s,".sam"))

    ## process remaining reads
    tmp<-readEval(sam, v)

    if(s == "T"){
      vcf$TreadSupport[which(vcf$ID == v)]<-paste0(sum(tmp == "WT read"),":",
                                                   sum(tmp == "Read supports only 1st position"),":",
                                                   sum(tmp == "Read supports only 2nd position"),":",
                                                   sum(tmp == "Mut read"))
    } else {
      if(s == "N"){
        vcf$NreadSupport[which(vcf$ID == v)]<-paste0(sum(tmp == "WT read"),":",
                                                     sum(tmp == "Read supports only 1st position"),":",
                                                     sum(tmp == "Read supports only 2nd position"),":",
                                                     sum(tmp == "Mut read"))
      }
    }
  }
  ## progress report
  print(paste0("Processed variant ",v," (",which(vcf$ID == v),"/",nrow(vcf),")"))
}
return(vcf)
}
# debug(tt)
vcf2 <- tt()
#use system2

# debug(DBSconclusion)
## Final conclusion on DBS evaluation
## I am allowing for a single read of either 1st or 2nd position supporting, to allow for a small number of sequencing errors
vcf<-DBSconclusion(vcf2)

write.table(vcf,file=paste0(smp,"_PCAWG_DBS_evaluated.vcf"),sep="\t",quote=F,row.names = F)

old <- data.table::fread("../MNVverify/data-raw/MELA-0161_PCAWG_DBS_evaluated.vcf")
cat("regression result", all.equal(old, data.table::data.table(vcf)))

## cleanup
tmp.files <- list.files(path = "./temp", full.names = TRUE)
if (unlink(tmp.files) == 1) message("unknown problem unlinking files in in ./temp")
print("Finished")
