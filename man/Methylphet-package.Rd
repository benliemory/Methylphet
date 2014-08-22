\name{Methylphet-package}
\alias{Methylphet-package}
\alias{Methylphet}
\docType{package}
\title{
Base-resolution methylation patterns accurately predict transcription factor bindings in vivo  

}
\description{
Methylphet adopt a novel supervised learning approach to predict TF-DNA interaction using DNA methylation profiles from whole-genome bisulfite sequencing data by exploiting the connection between DNA methylation level and TF binding. Two beta-binomial models are devised to characterize DNA methylation data around TF binding sites and the background. Along with other static genomic features, Methylphet adopt a supervised learning framework to predict TF-DNA interaction. When all methylation profile are taken together and combined with features at the sequence level, Methylphet can accurately predict TF binding and performs favorably when compared against competing methods.

}
\details{
\tabular{ll}{
Package: \tab Methylphet\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-08-21\cr
License: GPL-2
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
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
