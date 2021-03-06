% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VerifyDBSVcf.R
\name{VerifyDBSVcf}
\alias{VerifyDBSVcf}
\title{Given a VCF and the locations of a tumor and normal BAM files check whether analysis of individual reads supports DBSs calls.}
\usage{
VerifyDBSVcf(
  vcf,
  Nbam.name,
  Tbam.name,
  N.slice.dir,
  T.slice.dir,
  padding = 10,
  verbose = 0
)
}
\arguments{
\item{vcf}{An in-memory representation of a "variant call file", VCF, as a \code{data.frame}.}

\item{Nbam.name}{The name of the BAM file for the normal sample.}

\item{Tbam.name}{The name of the BAM file for the tumor sample.}

\item{N.slice.dir}{The directory containing the slices of the normal BAM.}

\item{T.slice.dir}{The directory containing the slices of the tumor BAM. Must
be different than \code{N.slice.dir}.}

\item{padding}{The number of base pairs on either side of the first position
of the DBS to include
in the slices.}

\item{verbose}{If > 0 print a message when starting the number of slices
generated every \code{verbose} slices.}
}
\value{
An in-memory VCF based on the input \code{vcf} with the additional columns
as described in \code{\link{Read_DBS_VCF_and_BAMs_to_verify_DBSs}}.
}
\description{
Given a VCF and the locations of a tumor and normal BAM files check whether analysis of individual reads supports DBSs calls.
}
\keyword{internal}
