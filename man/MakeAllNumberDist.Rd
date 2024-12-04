\name{MakeAllNumberDist}
\alias{MakeAllNumberDist}
\title{
Generates distributions for the counts of unique alleles in samples with different size g  (number of haploid genotypes).
}
\description{
By sampling with replacement from the allele frequency table for the reference population, with or without Fis correction, this function generates sampling distributions for the counts of unique alleles for samples of  different size g (number of haploid genotypes).
}
\usage{
MakeAllNumberDist(genotypes.df,ngametos=80,All.nchar=3,Fis=FALSE,Fis.treshold=0.1,
npermutations=5000) 
}
\arguments{
  \item{genotypes.df}{
The input file with the reference population microsatellite allele codes. This is a data.frame with one column per locus one row per individual with column names containing the locus names and without individual information (i.e. the first column is the information for the first loci). Use the same number of digits to code each allele. For example, 120126 is a possible value for a diploid coded with three digits per allele. At this stage, the reference data needs to be in diploid format. Missing data needs to be coded with 0, 00, 0000, or 000000..
}
 \item{ngametos}{
The method will generate sampling distributions for the counts of unique alleles for samples of  different size g (number of haploid genotypes). Distributions are sampled from the reference allele frequency table. There are distributions for each combination of locus and number of genotypes from 1 to \emph{ngametos}.
}
\item{All.nchar}{
How many digits are used in the reference population file to code each allele (can only be 1, 2, or 3).
}
\item{Fis}{
A logical indicating if Fis corrections should be applied when generating the sampling distributions. If yes, Genepop is used to calculate Qintra, the observed frequencies of identical pairs of genes in the reference population. Alleles are then sampled sequential with probability Qintra to sample the previously sampled allele and 1-Qintra to sample any other allele according to their allelic frequencies. Defaults to FALSE.
}
\item{Fis.treshold}{
When applying Fis correction, it is recommended to only correct for Fis values above a certain threshold. The default is 0.1.
}
 \item{npermutations}{
How many samples should each sampling distribution contain; defaults to 5000.
}
}
\details{
The probability \emph{P(OCUA | g)} that a given count of unique microstallite alleles \emph{OCUA}  is amplified at locus \emph{i} from an eDNA sample was produced by \emph{g} genotypes is estimated using a numerical maximum likelihod approach. The relative frequency of \emph{OCUA} in different sampling distributions of expected counts of unique allele (ECUA) is calculated for samples of size \emph{g}, from 1 to \emph{G}.  The multilocus \emph{P(OCUA | g)} is calculated by summing the log probabilities for all loci for each \emph{g} value. The maximum likelihood estimate is simply the maximum multilocus probability across the values of \emph{g}. The sampling distributions for each locus and \emph{g}  combination  are created by sampling with replacement samples of size \emph{g} alleles according to their frequencies in the reference population, with or without Fis correction,  and counting the number of unique alleles \emph{CUA} in each sample.
}
\value{
A list with two components. The first is a three-dimensional array with dimensions ngametos, number of loci (in the reference population allele frequency table), and npermutations, containing values for counts of unique alleles sampled. The second is a table with the allele frequencies in the reference population. This object will be required by the function \code{\link{G.Estimate}} to estimate the most likely number of genotypes present in the eDNA sample.
}
\references{
Liggan L, Alberto F (in prep)...
}
\author{
F. Alberto and L. Liggan.
}
\seealso{
\code{\link{G.Estimate}}
}
\keyword{ IO }
