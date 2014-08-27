\name{OtherGenomicFeatures.H1.chr10}
\alias{OtherGenomicFeatures.H1.chr10}
\docType{data}
\title{
Other Genomic Features on H1-hESC cell line(chromosome 10)
}
\description{
Other Genomic Features include distance between the center of candidate TFBS and the nearest TSS, distance between the center of candidate TFBS and the nearest center of other candidate TFBS, average conservation score of all bases on candidate TFBS, whether the candidate TFBS is overlapping with a repetitive region, whether the candidate TFBS is overlapping with an Intron, whether the candidate TFBS is overlapping with an Exon, whether the candidate TFBS is overlapping with an TSS, whether the candidate TFBS is overlapping with an TES, whether the candidate TFBS is overlapping with an CpG.

}
\usage{data(OtherGenomicFeatures.H1.chr10)}
\format{
  The format is:
 num [1:7233, 1:9] 0.3998 0.0767 0.0532 0.0173 0.0103 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:9] "motif.TSSdist" "motif.phastCon.score" "is.tss" "is.tes" ...
}
\examples{
data(OtherGenomicFeatures.H1.chr10)
## maybe str(OtherGenomicFeatures.H1.chr10) ; plot(OtherGenomicFeatures.H1.chr10) ...
}
\keyword{datasets}
