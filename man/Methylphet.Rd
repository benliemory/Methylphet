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
Candidate site location information and motif score for train data. This input is required.
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
Other Genomic features for train data. Genomic features include distance between the center of candidate TFBS and the nearest TSS, distance between the center of candidate TFBS and the nearest center of other candidate TFBS, average conservation score of all bases on candidate TFBS, etc.
}
  \item{testdata.mat}{
Candidate site location information and motif score for test data
}
  \item{testdata.methyl1}{
First set of methylation information for test data. Summarized methylation information +/- 10bins(30 bps/bin) around candidate sites. It could be CpG, CpH or 5hmc data. This input is required.
}
  \item{testdata.methyl2}{
Second set of methylation information for test data. Summarized methylation information +/- 10bins(30 bps/bin) around candidate sites. It could be CpG, CpH or 5hmc data. This input is optional.
}
  \item{OtherGenomicFeatures.test}{
Other Genomic features for test data. Genomic features include distance between the center of candidate TFBS and the nearest TSS, distance between the center of candidate TFBS and the nearest center of other candidate TFBS, average conservation score of all bases on candidate TFBS, etc.
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
Tianlei Xu, Ben Li, Hao Wu, Zhaohui Qin
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

### Load the package
library("Methylphet")


### Load all the data needed 
data(list = data(package="Methylphet")$results[,3])

### Using CpG methylation information only to predict TFBS when 0/1 golden standard is provided.
OCT4_CpG = Methylphet(traindata.mat = mESdata.motif.chr10,traindata.methyl1=mESdata.CpG.chr10, 
                      goldstandard=goldstandard.chr10, OtherGenomicFeatures.train=OtherGenomicFeatures.mES.chr10,
                      testdata.mat =H1data.motif.chr10, testdata.methyl1=H1data.CpG.chr10,
                      OtherGenomicFeatures.test=OtherGenomicFeatures.H1.chr10)

### Using both CpG and CpH methylation information to predict TFBS when location for ChIP-seq peaks are provided.
OCT4_CpG_CpH = Methylphet(traindata.mat = mESdata.motif.chr10,
                       traindata.methyl1=mESdata.CpG.chr10,traindata.methyl2=mESdata.CpH.chr10,
                          ChIPseqPeaks = peak.GR.ES.chr10, OtherGenomicFeatures.train=OtherGenomicFeatures.mES.chr10,
                      testdata.mat =H1data.motif.chr10, 
                          testdata.methyl1=H1data.CpG.chr10,testdata.methyl2=H1data.CpH.chr10,
                      OtherGenomicFeatures.test=OtherGenomicFeatures.H1.chr10)



}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
