\name{hmmbuild}
\alias{hmmbuild}
\title{
Create Profile Hidden Markov Model of given aligned sequences
}
\description{
This function creates profile hidden markov model of the given aligned sequences using HMMER algorithm.[1]
}
\usage{
hmmbuild(file_seq, file_out,pathHMMER="")
}
\arguments{
  \item{file_seq}{
Multiple sequence aligned fasta file
}
  \item{file_out}{
Output hidden markov model file
}
\item{pathHMMER}{
  Path where HMMER software is installed. Note: Windows user must setup cygwin to use this feature and set path to HMMER binaries ( ~hmmer/binaries/)
}
}
\value{
Create Profile Hidden Markov Model in local directory
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
\code{\link{nhmmer}}
}
\examples{
msa = system.file("extdata", "msa.fasta", package = "MethTargetedNGS")
if (file.exists("/usr/bin/hmmbuild"))
  hmmbuild(file_seq=msa,file_out="hmm",pathHMMER = "/usr/bin")
}
\keyword{ ProfileHMM }
\keyword{ hmmbuild }
