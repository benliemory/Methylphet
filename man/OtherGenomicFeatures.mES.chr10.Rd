\name{OtherGenomicFeatures.mES.chr10}
\alias{OtherGenomicFeatures.mES.chr10}
\docType{data}
\title{
Other Genomic Features
}
\description{
Other Genomic Features include distance between the center of candidate TFBS and the nearest TSS, distance between the center of candidate TFBS and the nearest center of other candidate TFBS, average conservation score of all bases on candidate TFBS, whether the candidate TFBS is overlapping with a repetitive region, whether the candidate TFBS is overlapping with an Intron, whether the candidate TFBS is overlapping with an Exon, whether the candidate TFBS is overlapping with an TSS, whether the candidate TFBS is overlapping with an TES, whether the candidate TFBS is overlapping with an CpG.
}
\usage{data(OtherGenomicFeatures.mES.chr10)}
\format{
  The format is:
 num [1:10973, 1:9] 0.00773 0.00873 0.0098 0.02461 0.03432 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:9] "motif.TSSdist" "motif.phastCon.score" "is.tss" "is.tes" ...
}

\examples{
data(OtherGenomicFeatures.mES.chr10)
}
\keyword{datasets}
