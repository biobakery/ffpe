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
Sample metadata includes input RNA concentration, cell type, and
replicate #.

}

\usage{data(lumibatch.GSE17565)}
\format{
  
  Formal class 'LumiBatch' [package "lumi"]
  
}

\details{
  
  Downloaded from the GEO URL:
  
  http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE17565 .

  Raw data were obtained from the supplemental file
  GSE17565_nonorm_nobkgd.txt.gz, and inserted into a lumibatch object
  using the lumi:lumiR command.
  
}

\source{

  April C, Klotzle B, Royce T, Wickham-Garcia E et al. Whole-genome gene
  expression profiling of formalin-fixed, paraffin-embedded tissue
  samples. PLoS One 2009 Dec 3;4(12):e8162. PMID: 19997620
  
  Data from accession ID GSE17565 of the Gene Expression Omnibus
  
}

\examples{

library(ffpe)
data(lumibatch.GSE17565)
pData(lumibatch.GSE17565)
sortedIqrPlot(lumibatch.GSE17565,main="GSE17565")

}

\keyword{datasets}
