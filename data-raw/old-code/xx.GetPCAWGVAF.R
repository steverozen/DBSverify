xx.GetPCAWGVAF<-function(vcf){
  tmp<-vcf[[8]]
  tmp<-stringr::str_split_fixed(tmp,"VAF=",2)[,2]
  tmp<-stringr::str_split_fixed(tmp,";",2)[,1]
  tmp<-as.numeric(as.character(tmp))
  return(cbind(vcf, VAF=tmp, read.depth = rep(0,length(tmp))))
}
