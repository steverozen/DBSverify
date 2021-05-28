ReadVCFAndBAMsAndProcess <- function(vcf.file, norm) {
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

  # write.table(vcf.list$DBS[[1]],file=paste0(smp,"_PCAWG_DBS.vcf"),sep="\t",row.names = F,quote=F)

  ## read DBS ICAMS vcf
  # vcf<-read.csv(paste0(smp,"_PCAWG_DBS.vcf"),sep="\t",as.is=T)

  # vcf <- vcf.list$DBS[[1]]

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
  tt <- function(vcf) {
    for(v in 1:nrow(vcf)){

      ## first create bam slices
      ## add a slight window around the variant, variants that don't overlap the (entire) DBS will be ignored anyway
      # bamCoord <- paste0(vcf[v,"CHROM"],":",vcf[v,"POS"]-10,"-",vcf[v,"POS"]+10)
      # system(paste0("samtools view -h ",Nbam," ",bamCoord," > temp/",smp,"_",bamCoord,"_N.sam"),wait=T)
      # system(paste0("samtools view -h ",Tbam," ",bamCoord," > temp/",smp,"_",bamCoord,"_T.sam"),wait=T)
      CHROM = vcf[v, "CHROM"]
      POS   = vcf[v, "POS"]
      REF   = vcf[v, "REF"]
      ALT   = vcf[v, "ALT"]
      vcf[v, "TreadSupport"] <-
        GetReadSupport(BAM.name = Tbam, CHROM = CHROM, POS = POS, REF = REF, ALT = ALT)
      vcf[v, "NreadSupport"] <-
        GetReadSupport(BAM.name = Nbam, CHROM = CHROM, POS = POS, REF = REF, ALT = ALT)

      ## progress report
      # print(paste0("Processed variant ",v," (",which(vcf$ID == v),"/",nrow(vcf),")"))
    }
    return(vcf)
  }

  vcf2 <- tt(vcf)
  #use system2

  # debug(DBSconclusion)
  ## Final conclusion on DBS evaluation
  ## I am allowing for a single read of either 1st or 2nd position supporting, to allow for a small number of sequencing errors
  vcf3 <-DBSconclusion(vcf2)

  write.table(vcf3, file=paste0(smp,"_PCAWG_DBS_evaluated.vcf"), sep="\t", quote=F, row.names = F)

  old <- data.table::fread("../MNVverify/data-raw/MELA-0161_PCAWG_DBS_evaluated.vcf")
  cat("regression result", all.equal(old, data.table::data.table(vcf3)))


  print("Finished")




}
