\name{mESdata.motif.chr10}
\alias{mESdata.motif.chr10}
\docType{data}
\title{
OCT4 mm9 motif(chr10) 
}
\description{
candidate location and motif score.
}
\usage{data(mESdata.motif.chr10)}
\format{
  A data frame with 10973 observations on the following 6 variables.
  \describe{
    \item{\code{seqnames}}{a factor with levels \code{chr1} \code{chr2} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chr10} \code{chr11} \code{chr12} \code{chr13} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chrX} \code{chrY}}
    \item{\code{start}}{a numeric vector}
    \item{\code{end}}{a numeric vector}
    \item{\code{width}}{a numeric vector}
    \item{\code{strand}}{a factor with levels \code{+} \code{-} \code{*}}
    \item{\code{motif}}{a numeric vector}
  }
}

\examples{
data(mESdata.motif.chr10)
## maybe str(mESdata.motif.chr10) ; plot(mESdata.motif.chr10) ...
}
\keyword{datasets}
