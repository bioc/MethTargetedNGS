\name{bconv}
\alias{bconv}
\title{
Convert non-bisulfite sequences to bisulfite sequences.
}
\description{
Bisulfite sequences are the bisulfite treated DNA sequences where all cytosines except cytosine from CpG sites are converted to thymie. This technique is used to determine pattern of methylation. This function convert all cytosine except cytosines from CpG sites to thymine.
}
\usage{
bconv(fasta_file, out_file = "output.fasta")
}
\arguments{
  \item{fasta_file}{
   : Input fasta file for conversion
}
  \item{out_file}{
    : String value naming an output file. Default is output.fasta   
}
}
\value{
Fasta File
}
\author{
Muhammad Ahmer Jamil, Prof. Holger Frohlich, Priv.-Doz. Dr. Osman El-Maarri

Maintainer: Muhammad Ahmer Jamil \email{engr.ahmerjamil@gmail.com}
}
\examples{
input =  system.file("extdata", "bisulfite.fasta", package = "MethTargetedNGS")
bconv(fasta_file = input, out_file = "output.fasta")
}
\keyword{ Bisulfite Conversion}
\keyword{ Methylation }