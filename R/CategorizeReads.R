#' Categorizes the reads in one sam file.
#'
#' @param sam An in-memory version of a sam file as data.frame.
#'
#' @param POS The position in the sam at which the DBS starts.
#'
#' @param REF The reference allele.
#'
#' @param ALT The putative alternative allele.
#'
#' @return a named character vector, the names are the read names in the sam file.
#'  Each element is one of
#'  "Read supports only 1st position"
#'  "Read supports only 2nd position",
#   "Read does not overlap complete DBS",
#'  "WT read"
#'  "Mut read"
#'
#'  @keywords internal
#'

CategorizeReads <-function(sam, POS, REF, ALT){
  POS <- as.numeric(POS) # In case it is a character string
  tmp<-NULL
  if (nrow(sam) == 0) { stop("Do not call CategorizeReads when nrow(sam) == 0")}
  for(r in 1:nrow(sam)){ # every row is a read
    ## check whether the entire DBS is covered by the read
    readRange<-c(as.numeric(sam$POS[r]),as.numeric(sam$POS[r])+nchar(sam$SEQ[r])-1)
    if(!(POS >= min(readRange) & POS+1 <= max(readRange))){
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

      readSeq<-substr(sam[r,"SEQ"], POS-as.numeric(sam[r,"POS"])+1,
                      POS-as.numeric(sam[r,"POS"])+1+1)

      ## check whether the sequence of the read matches the DBS
      pos1eval<-substr(readSeq,1,1) == substr(REF,1,1)
      pos2eval<-substr(readSeq,2,2) == substr(REF,2,2)
      if(pos1eval && pos2eval){
        tmp<-c(tmp,"WT read")
        names(tmp)[length(tmp)]<-sam$QNAME[r]
      } else {
        pos1eval<-substr(readSeq,1,1) == substr(ALT, 1,1)
        pos2eval<-substr(readSeq,2,2) == substr(ALT ,2,2)
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
