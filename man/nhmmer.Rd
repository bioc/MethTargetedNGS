\name{nhmmer}
\alias{nhmmer}
\title{
Calculate likelihood of the given profile hidden markov model against group of sequences
}
\description{
This function calculates likelihood score of given pool of sequences against given profile hidden markov model using HMMER algorithm.[1]
}
\usage{
nhmmer(file_hmm, file_seq, pathHMMER="")
}
\arguments{
  \item{file_hmm}{
  HMM file from hmmbuild function
}
  \item{file_seq}{
  Sequence fasta file for calculating likelihood
}
  \item{pathHMMER}{
   Path where HMMER software is installed. Note: Windows user must setup cygwin to use this feature and set path to HMMER binaries ( ~hmmer/binaries/)
  }
}
\value{
Matrix containing likelihood scores
}
\references{
[1]Finn, Robert D., Jody Clements, and Sean R. Eddy. "HMMER web server: interactive sequence similarity searching." Nucleic acids research (2011): gkr367.
}
\author{
Muhammad Ahmer Jamil, Prof. Holger Frohlich, Priv.-Doz. Dr. Osman El-Maarri

Maintainer: Muhammad Ahmer Jamil \email{engr.ahmerjamil@gmail.com}
}
\note{
Require HMMER software

Windows User: Please download HMMER from http://hmmer.janelia.org/

Setup cygwin from http://www.cygwin.com

Linux/Mac User: Download binaries or compile HMMER from http://hmmer.janelia.org/
}
\seealso{
\code{\link{hmmbuild}}
}
\examples{
msa = system.file("extdata", "msa.fasta", package = "MethTargetedNGS")
tumor = system.file("extdata", "Tumor.fasta", package = "MethTargetedNGS")
if (file.exists("/usr/bin/hmmbuild"))
{hmmbuild(file_seq=msa,file_out="hmm",pathHMMER = "/usr/bin")
res <- nhmmer("hmm",tumor,pathHMMER = "/usr/bin")
res}
}
\keyword{ HMMER }
\keyword{ ProfileHMM }
\keyword{ Methylation }
