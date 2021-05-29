SummarizeReadSupport <- function(vcf, Nbam.name, Tbam.name) {
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

    vcf[v, "NreadSupport"] <-
      GetReadSupport(BAM.name = Nbam.name, CHROM = CHROM, POS = POS, REF = REF, ALT = ALT)

    vcf[v, "TreadSupport"] <-
      GetReadSupport(BAM.name = Tbam.name, CHROM = CHROM, POS = POS, REF = REF, ALT = ALT)
  }
  return(vcf)
}
