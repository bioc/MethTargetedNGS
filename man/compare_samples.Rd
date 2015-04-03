\name{compare_samples}
\alias{compare_samples}
\title{
 Complete Methylation Analysis of Next Generation Sequencing Data
}
\description{
This function perform complete methylation analysis of the data. 

1. Visualize methylation pattern

2. Calculate methylation average

3. Calculate methylation entropy

4. Perform fisher exact test on the samples to identify significant CpG sites.
}
\usage{
compare_samples(healthy, tumor)
}
\arguments{
  \item{healthy}{
 : Output Matrix from \code{methAlign}
}
  \item{tumor}{
: Output Matrix from \code{methAlign}
}
}
\value{
Generate a plot of Methylation Average, Methylation Entropy, Fisher Exact Test and Log Odd Ratio
}
\author{
Muhammad Ahmer Jamil, Prof. Holger Frohlich, Priv.-Doz. Dr. Osman El-Maarri

Maintainer: Muhammad Ahmer Jamil \email{engr.ahmerjamil@gmail.com}
}
\note{
This function needs time to process depending on the number of sequences in fasta file
}

\seealso{
\code{\link{methAlign}},
\code{\link{methAvg}},
\code{\link{methEntropy}},
\code{\link{odd_ratio}},
\code{\link{fishertest_cpg}},
}
\examples{
healthy = system.file("extdata", "Healthy.fasta", package = "MethTargetedNGS")
tumor = system.file("extdata", "Tumor.fasta", package = "MethTargetedNGS")
reference =  system.file("extdata", "Reference.fasta", package = "MethTargetedNGS")

healthy = methAlign(healthy,reference)
tumor = methAlign(tumor,reference)
compare_samples(healthy,tumor)
}
\keyword{ Methylation }
\keyword{ Average }
\keyword{ Entropy }
\keyword{ Odd Ratio }
\keyword{ Significant CpG sites}
\keyword{ Fisher Test }