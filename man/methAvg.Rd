\name{methAvg}
\alias{methAvg}
\title{
Calculate Methylation Average of given methylation matrix
}
\description{
Methylation average of a CpG site is the percentage of unmethylated cytosine or methylated cytosine in a particular CpG site. The methylation average of a particular CpG site was calculated by number of cytosine divided by sum of total number of methylated and unmethylated cytosine at particular CpG site in a group of reads.

            average = NC/(NC + NT)

}
\usage{
methAvg(Sample, plot = FALSE)
}

\arguments{
  \item{Sample}{
  Matrix from methAlign. Also matrix where columns represents Cytosine of CpG sites and rows represents sequences.

}
  \item{plot}{
  Boolean. TRUE if need a plot after calculation. Default FALSE

}
}
\value{
 Vector containing average methylation of given methylation matrix. Length of the vector represents the number of CpG sites in methylation matrix.
}
\author{
Muhammad Ahmer Jamil, Prof. Holger Frohlich, Priv.-Doz. Dr. Osman El-Maarri

Maintainer: Muhammad Ahmer Jamil \email{engr.ahmerjamil@gmail.com}
}
\seealso{
\code{\link{methAlign}},
\code{\link{compare_samples}}
}
\examples{
healthy = system.file("extdata", "Healthy.fasta", package = "MethTargetedNGS")
reference = system.file("extdata", "Reference.fasta", package = "MethTargetedNGS")
methP <- methAlign(healthy,reference)
avgMeth <- methAvg(methP,plot=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Average}
\keyword{Methylation}
% __ONLY ONE__ keyword per line
