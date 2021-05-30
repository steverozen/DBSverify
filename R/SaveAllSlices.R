SaveAllSlices <- function(vcf, bam.name, where.to.put.slices) {
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

    file.name <-
      file.path(where.to.put.slices,
                paste("...part of BAM name...", CHROM, POS, REF, ALT, sep = "-"))
    # Create saved BAM slices
    CreateBAMSliceFileSamtools(BAM.name, CHROM, POS, padding = 10, save.file.path = file.name)
  }

}
