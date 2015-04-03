\name{odd_ratio}
\alias{odd_ratio}
\title{
Calculate log odd ratio of the given samples (healthy/tumor)
}
\description{
Log Odd ratio defines the hypomethylation and hypermethylation of a sample in comparison to the other sample.
}
\usage{
odd_ratio(SampA, SampB, plot = TRUE, main = "Log Odd Ratio")
}
\arguments{
  \item{SampA}{
Matrix from methAlign. Also matrix where columns represents Cytosine of CpG sites and rows represents sequences.
}
  \item{SampB}{
Matrix from methAlign. Also matrix where columns represents Cytosine of CpG sites and rows represents sequences.
}
  \item{plot}{
Boolean. TRUE if need a plot after calculation. Default TRUE
}
  \item{main}{
Title of the plot
}
}
\value{
Vector containing log odd ratios.
}
\author{
Muhammad Ahmer Jamil, Prof. Holger Frohlich, Priv.-Doz. Dr. Osman El-Maarri

Maintainer: Muhammad Ahmer Jamil \email{engr.ahmerjamil@gmail.com}
}
\seealso{
\code{\link{methAlign}}
}
\examples{
healthy = system.file("extdata", "Healthy.fasta", package = "MethTargetedNGS")
tumor = system.file("extdata", "Tumor.fasta", package = "MethTargetedNGS")
reference =  system.file("extdata", "Reference.fasta", package = "MethTargetedNGS")

healthy = methAlign(healthy,reference)
tumor = methAlign(tumor,reference)
odd_ratio(healthy,tumor)
}
\keyword{ OddRatio }
\keyword{ Methylation }