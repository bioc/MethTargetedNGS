\name{methHeatmap}
\alias{methHeatmap}
\title{
Generate Heatmap of the given methylation data.
}
\description{
Heatmaps are the way of visualizing methylation statuses of a sample. This function allows user to visualize methylation statuses at each CpG site for every sequence available in pool.
}
\usage{
methHeatmap(Sample, yl = "", plot = TRUE, title = "")
}

\arguments{
  \item{Sample}{
Matrix from methAlign. Also matrix where columns represents Cytosine of CpG sites and rows represents sequences.
}
  \item{yl}{
Ylabel for heatmap
}
  \item{plot}{
Boolean. If plot == FALSE, function will return a matrix of 1s and 0s. If plot == TRUE, function will create a heatmap as well as return a matrix of 1s and 0s
}
  \item{title}{
Title of the heatmap
}
}
\value{
Heatmap
}
\author{
Ahmer Jamil \email{engr.ahmerjamil@gmail.com}
}
\seealso{
\code{\link{methAlign}}
}
\examples{
healthy = system.file("extdata", "Healthy.fasta", package = "MethTargetedNGS")
reference = system.file("extdata", "Reference.fasta", package = "MethTargetedNGS")
healthy = methAlign(healthy,reference)
hHeatmap = methHeatmap(healthy,plot=TRUE)
}
\keyword{ Heatmap }
\keyword{ Methylation }