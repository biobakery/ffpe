\name{lumibatch.GSE17565}
\alias{lumibatch.GSE17565}
\docType{data}
\title{
  Illumina DASL expression data from FFPE tissues.
}
\description{
This lumibatch object contains raw expression data for 32 FFPE samples
of Burkitts Lymphoma and Breast Adenocarcinoma, with dilution series and
technical duplicates.  From the original study by April et al. (2009).
Sample metadata includes input RNA concentration, cell type, and replicate #.
}
\usage{data(lumibatch.GSE17565)}
\format{
  The format is:
Formal class 'LumiBatch' [package "lumi"] with 10 slots
  ..@ history          :'data.frame':	3 obs. of  4 variables:
  .. ..$ submitted  : chr [1:3] "2011-10-21 19:26:35" "2011-10-21 19:26:41" "2011-10-21 19:46:14"
  .. ..$ finished   : chr [1:3] "2011-10-21 19:26:41" "2011-10-21 19:26:42" "2011-10-21 19:46:14"
  .. ..$ command    : chr [1:3] "lumiR(\"GSE17565_nonorm_nobkgd.txt\", inputAnnotation = FALSE)" "lumiQ(x.lumi = x.lumi, detectionTh = detectionTh, verbose = verbose)" "Subsetting 32 samples."
  .. ..$ lumiVersion: chr [1:3] "2.2.1" "2.2.1" "2.2.1"
  ..@ controlData      :'data.frame':	0 obs. of  0 variables
  ..@ QC               :List of 3
  .. ..$ BeadStudioSummary:'data.frame':	0 obs. of  0 variables
  .. ..$ sampleSummary    : num [1:4, 1:32] 9.264 3.151 0.536 225.7 9.182 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "mean" "standard deviation" "detection rate(0.01)" "distance to sample mean"
  .. .. .. ..$ : chr [1:32] "GSM437726" "GSM437727" "GSM437728" "GSM437729" ...
  .. ..$ history          :'data.frame':	2 obs. of  4 variables:
  .. .. ..$ submitted  : chr [1:2] "2011-10-21 19:26:35" "2011-10-21 19:26:41"
  .. .. ..$ finished   : chr [1:2] "2011-10-21 19:26:41" "2011-10-21 19:26:42"
  .. .. ..$ command    : chr [1:2] "lumiR(\"GSE17565_nonorm_nobkgd.txt\", inputAnnotation = FALSE)" "lumiQ(x.lumi = x.lumi, detectionTh = detectionTh, verbose = verbose)"
  .. .. ..$ lumiVersion: chr [1:2] "2.2.1" "2.2.1"
  ..@ assayData        :<environment: 0x221aaf8> 
  ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	4 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr [1:4] NA NA NA NA
  .. .. ..@ data             :'data.frame':	32 obs. of  4 variables:
  .. .. .. ..$ title      : chr [1:32] "Raji 100 ng rep1" "Raji 100 ng rep2" "Raji 50 ng rep1" "Raji 50 ng rep2" ...
  .. .. .. ..$ inputRNA.ng: num [1:32] 100 100 50 50 10 10 5 5 1 1 ...
  .. .. .. ..$ cell.type  : chr [1:32] "Burkitts Lymphoma" "Burkitts Lymphoma" "Burkitts Lymphoma" "Burkitts Lymphoma" ...
  .. .. .. ..$ replicate  : chr [1:32] "1" "2" "1" "2" ...
  .. .. ..@ dimLabels        : chr [1:2] "rowNames" "columnNames"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	1 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr "The Illumina microarray identifier"
  .. .. ..@ data             :'data.frame':	24526 obs. of  1 variable:
  .. .. .. ..$ ID_REF: chr [1:24526] "ILMN_1762337" "ILMN_2055271" "ILMN_2383229" "ILMN_1806310" ...
  .. .. ..@ dimLabels        : chr [1:2] "featureNames" "featureColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots
  .. .. ..@ name             : chr ""
  .. .. ..@ lab              : chr ""
  .. .. ..@ contact          : chr ""
  .. .. ..@ title            : chr ""
  .. .. ..@ abstract         : chr ""
  .. .. ..@ url              : chr ""
  .. .. ..@ pubMedIds        : chr ""
  .. .. ..@ samples          : list()
  .. .. ..@ hybridizations   : list()
  .. .. ..@ normControls     : list()
  .. .. ..@ preprocessing    : list()
  .. .. ..@ other            :List of 1
  .. .. .. ..$ Data File Information: chr "# GEOarchive matrix file."
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 0 0
  ..@ annotation       : chr(0) 
  ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr(0) 
  .. .. ..@ data             :'data.frame':	32 obs. of  0 variables
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. ..@ .Data:List of 4
  .. .. .. ..$ : int [1:3] 2 12 1
  .. .. .. ..$ : int [1:3] 2 10 0
  .. .. .. ..$ : int [1:3] 1 3 0
  .. .. .. ..$ : int [1:3] 1 0 0
}
\details{
  
  Downloaded from the GEO URL:
  
  http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE17565 .

  Raw data were obtained from the supplemental file
  GSE17565_nonorm_nobkgd.txt.gz, and inserted into a lumibatch object
  using the lumi:lumiR command.
}

\source{

  April C, Klotzle B, Royce T, Wickham-Garcia E et al. Whole-genome gene expression profiling of formalin-fixed, paraffin-embedded tissue samples. PLoS One 2009 Dec 3;4(12):e8162. PMID: 19997620

}
\examples{

library(ffpe)
data(lumibatch.GSE17565)
pData(lumibatch.GSE17565)
sortedIqrPlot(lumibatch.GSE17565,main="GSE17565")

}
\keyword{datasets}
