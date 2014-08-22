\name{Methylphet}
\alias{Methylphet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methylphet
}
\description{
Methylphet main function to predict transcription factor bindings in vivo 
}
\usage{
Methylphet(traindata.mat = NA, traindata.methyl1 = NA, traindata.methyl2 = NA, goldstandard = NA, ChIPseqPeaks = NA, OtherGenomicFeatures.train = NA, testdata.mat = NA, testdata.methyl1 = NA, testdata.methyl2 = NA, OtherGenomicFeatures.test = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traindata.mat}{
candidate site location and motif score for train data
}
  \item{traindata.methyl1}{
First set of methylation information for train data. Summarized methylation information +/- 10bins(30 bps/bin) around candidate sites. It could be CpG, CpH or 5hmc data. This input is required.
}
  \item{traindata.methyl2}{
Second set of methylation data. Summarized methylation information +/- 10bins(30 bps/bin) around candidate sites. It could be CpG, CpH or 5hmc data. This input is optional.
}
  \item{goldstandard}{
0/1 results indicating whether the candidate sites are overlap with ChIP-seq peaks. Either goldstandard or ChIPseqPeaks are required.
}
  \item{ChIPseqPeaks}{
Location information for ChIP-seq peaks. Either goldstandard or ChIPseqPeaks are required.
}
  \item{OtherGenomicFeatures.train}{
Other Genomic features including conservation scores, etc. for train data.
}
  \item{testdata.mat}{
candidate site location and motif score for test data
}
  \item{testdata.methyl1}{
First set of methylation information for test data. Summarized methylation information +/- 10bins(30 bps/bin) around candidate sites. It could be CpG, CpH or 5hmc data. This input is required.
}
  \item{testdata.methyl2}{
Second set of methylation information for test data. Summarized methylation information +/- 10bins(30 bps/bin) around candidate sites. It could be CpG, CpH or 5hmc data. This input is optional.
}
  \item{OtherGenomicFeatures.test}{
Other Genomic features including conservation scores, etc.for test data.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
data()

OCT4_CpG = Methylphet(traindata.mat = mESdata.motif.chr10,traindata.methyl1=mESdata.CpG.chr10, 
                      goldstandard=goldstandard.chr10, OtherGenomicFeatures.train=OtherGenomicFeatures.mES.chr10,
                      testdata.mat =H1data.motif.chr10, testdata.methyl1=H1data.CpG.chr10,
                      OtherGenomicFeatures.test=OtherGenomicFeatures.H1.chr10)


OCT4_CpG_CpH = Methylphet(traindata.mat = mESdata.motif.chr10,
                       traindata.methyl1=mESdata.CpG.chr10,traindata.methyl2=mESdata.CpH.chr10,
                          ChIPseqPeaks = peak.GR.ES.chr10, OtherGenomicFeatures.train=OtherGenomicFeatures.mES.chr10,
                      testdata.mat =H1data.motif.chr10, 
                          testdata.methyl1=H1data.CpG.chr10,testdata.methyl2=H1data.CpH.chr10,
                      OtherGenomicFeatures.test=OtherGenomicFeatures.H1.chr10)



}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
