\name{methEntropy}
\alias{methEntropy}
\title{
Calculate Methylation Entropy
}
\description{
Entropy comparison between healthy and tumor samples can identify significant CpG sites which are contributing most in the tumor development either by hypomethylation or hypermethylation. Also such way can help in understanding the randomness in methylation status.
Sliding window of 4 was used to calculate the entropy in the sample, which can analyze 16 different pattern for entropy calculation.

}
\usage{
methEntropy(x)
}
\arguments{
  \item{x}{
   Matrix from methAlign. Also matrix where columns represents Cytosine of CpG sites and rows represents sequences
}
}
\value{
Matrix containing entropy for every sequence and group of 4 cpg sites.
}
\references{
Xie, H., Wang, M., de Andrade, A., Bonaldo, M.d.F., Galat, V., Arndt, K., Rajaram, V.,
Goldman, S., Tomita, T. and Soares, M.B. (2011) Genome-wide quantitative assessment of
variation in DNA methylation patterns. Nucleic Acids Research, 39, 4099-4108.
}
\author{
Muhammad Ahmer Jamil, Prof. Holger Frohlich, Priv.-Doz. Dr. Osman El-Maarri

Maintainer: Muhammad Ahmer Jamil \email{engr.ahmerjamil@gmail.com}
}
\note{
This function needs time to process depending on the number of rows in matrix
}
\seealso{
\code{\link{methAlign}}
}
\examples{
healthy = system.file("extdata", "Healthy.fasta", package = "MethTargetedNGS")
reference = system.file("extdata", "Reference.fasta", package = "MethTargetedNGS")
methP <- methAlign(healthy,reference)
entMeth <- methEntropy(methP)
plot(entMeth,type="l")
}
\keyword{ Entropy }
\keyword{ Methylation Entropy }