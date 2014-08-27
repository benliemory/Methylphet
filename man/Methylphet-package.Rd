\name{Methylphet-package}
\alias{Methylphet-package}
\alias{Methylphet}
\docType{package}
\title{
Base-resolution methylation patterns accurately predict transcription factor bindings in vivo  
}
\description{
Methylphet adopts a two-step method to predict transcription factor-DNA interaction using DNA methylation profiles from whole-genome bisulfite sequencing data by exploiting the connection between DNA methylation level and transcription factor binding. In the first step, beta-binomial models are devised to characterize DNA methylation data around TF binding sites and the background to estimate methylation scores. Along with other static genomic features, a random forest framework is adopted in the second step to predict transcription factor-DNA interaction. When all methylation profile are taken together and combined with features at the sequence level, Methylphet can accurately predict TF binding and performs favorably when compared against competing methods.

}
\details{
\tabular{ll}{
Package: \tab Methylphet\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-08-21\cr
License: \tab GPL-2
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Tianlei Xu, Ben Li, Hao Wu, Zhaohui Qin

Maintainer: Tianlei Xu<tianlei.xu@emory.edu>, Ben Li<ben.li@emory.edu>, Hao Wu<hao.wu@emory.edu>, Zhaohui Qin<zhaohui.qin@emory.edu>

}


\keyword{ package }

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
