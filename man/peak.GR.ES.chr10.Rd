\name{peak.GR.ES.chr10}
\alias{peak.GR.ES.chr10}
\docType{data}
\title{
OCT4 mESC ChIP-seq Peaks
}
\description{
Location information for ChIP-seq peaks. This can convert to golden standard which is required for the train data. 
}
\usage{data(peak.GR.ES.chr10)}
\format{
  The format is:
Formal class 'GRanges' [package "GenomicRanges"] with 6 slots
  ..@ seqnames       :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. ..@ values         : Factor w/ 22 levels "chr1","chr10",..: 2
  .. .. ..@ lengths        : int 521
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
  .. .. ..@ start          : int [1:521] 3134809 3557539 3796733 4181829 4488304 5731965 5970992 6012984 6056693 6091822 ...
  .. .. ..@ width          : int [1:521] 159 159 176 256 245 267 320 168 168 159 ...
  .. .. ..@ NAMES          : NULL
  .. .. ..@ elementType    : chr "integer"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ strand         :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 3
  .. .. ..@ lengths        : int 521
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ elementMetadata:Formal class 'DataFrame' [package "IRanges"] with 6 slots
  .. .. ..@ rownames       : NULL
  .. .. ..@ nrows          : int 521
  .. .. ..@ listData       : Named list()
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ seqinfo        :Formal class 'Seqinfo' [package "GenomicRanges"] with 4 slots
  .. .. ..@ seqnames   : chr [1:22] "chr1" "chr10" "chr11" "chr12" ...
  .. .. ..@ seqlengths : int [1:22] NA NA NA NA NA NA NA NA NA NA ...
  .. .. ..@ is_circular: logi [1:22] NA NA NA NA NA NA ...
  .. .. ..@ genome     : chr [1:22] NA NA NA NA ...
  ..@ metadata       : list()
}
\examples{
data(peak.GR.ES.chr10)
## maybe str(peak.GR.ES.chr10) ; plot(peak.GR.ES.chr10) ...
}
\keyword{datasets}
