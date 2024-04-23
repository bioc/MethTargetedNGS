\name{methAlign}
\alias{methAlign}
\title{
 Align sequences with the reference sequence using pairwiseAlignment function from the pwalign package
}
\description{
This function allow users to align pool of sequences to the reference sequence.
}
\usage{
methAlign(sequence_fasta, ref_seq, sub_mat = FALSE, align_type = "local")
}
\arguments{
  \item{sequence_fasta}{
String value naming an input fasta file. Single sequence or Multiple sequences in a single fasta file
}
  \item{ref_seq}{
 String value naming an input fasta file. Single reference sequence is requried. If multiple sequences were passed only first sequence will be considered as reference.
}
  \item{sub_mat}{
 Substitution matrix for the alignment.
}
  \item{align_type}{
 type of alignment. One of "global", "local", "overlap", "global-local", and "local-global" where 
 
 "global" = align whole strings with end gap penalties, 
 
 "local" = align string fragments, 
 
 "overlap" = align whole strings without end gap penalties, 
 
 "global-local" = align whole strings in pattern with consecutive subsequence of subject, 
 
 "local-global" = align consecutive subsequence of pattern with whole strings in subject.
 
 Default is "local"
}
}
\value{
 Methylation Matrix. Number of rows represents number of reads in sequence fasta file and number of columns represents number of CpG sites in reference fasta sequence. Only Cytosine of CpG site was observed in the table whether it is methylated or unmethylated.

}
\author{
Muhammad Ahmer Jamil, Prof. Holger Frohlich, Priv.-Doz. Dr. Osman El-Maarri

Maintainer: Muhammad Ahmer Jamil \email{engr.ahmerjamil@gmail.com}
}
\note{
This function need some time to process depending on the number of sequences in fasta file
}
\seealso{
\code{\link{compare_samples}}
}
\examples{
healthy = system.file("extdata", "Healthy.fasta", package = "MethTargetedNGS")
reference = system.file("extdata", "Reference.fasta", package = "MethTargetedNGS")
methAlign(healthy,reference)
}
\keyword{ Sequence Alignment }
\keyword{ Methylation Analysis }
