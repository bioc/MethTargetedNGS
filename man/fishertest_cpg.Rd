\name{fishertest_cpg}
\alias{fishertest_cpg}
\title{
Perform Fisher Exact Test on Methylation Data.
}
\description{
Fisher exact test is a test to calculate the statistical significance using contingency table. It was used to find the statistically significant differences in the methylation status of one particular CpG site between healthy and tumor sample. Contingency matrix was created for each CpG site. P-value was corrected for multiple testing using Benjamini-Hochberg method to calculate False Discovery Rate (FDR)
}
\usage{
fishertest_cpg(healthy, tumor, plot = TRUE, main = "Fisher Exact Test")
}
\arguments{
  \item{healthy}{
Matrix from methAlign. Also matrix where columns represents Cytosine of CpG sites and rows represents sequences.
}
  \item{tumor}{
Matrix from methAlign. Also matrix where columns represents Cytosine of CpG sites and rows represents sequences.
}
  \item{plot}{
Boolean. TRUE if need a plot after calculation. Default TRUE
}
  \item{main}{
Title of the plot. Default "Fisher Exact Test"
}
}
\value{
Vector containing p-values.
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
tumor = system.file("extdata", "Tumor.fasta", package = "MethTargetedNGS")
reference =  system.file("extdata", "Reference.fasta", package = "MethTargetedNGS")

healthy = methAlign(healthy,reference)
tumor = methAlign(tumor,reference)
fisherexacttest <- fishertest_cpg(healthy,tumor)
}
\keyword{ Significant CpG }
\keyword{ Methylation }
\keyword{ Fisher Test}
