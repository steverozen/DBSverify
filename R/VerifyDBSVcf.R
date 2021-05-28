VerifyDBSVcf <- function(vcf, Nbam, Tbam) {

  # debug(addID)
  vcf$ID <- IDVector(vcf)
  ## remove duplicates if present
  vcf<-vcf[!duplicated(vcf$ID), ]
  rownames(vcf)<-vcf$ID

  vcf$Format<-paste0("WtReads:pos1reads:pos2reads:MutReads")

  vcf2 <- SummarizeReadSupport(vcf = vcf, Nbam = Nbam, Tbam = Tbam)
  #use system2

  vcf3 <-DBSConclusion(vcf2)
  return(vcf3)
}
