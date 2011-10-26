\name{sortedIqrPlot}

\title{
  Minimal box plot of samples, which get ordered optionally by batch,
  then by IQR.
}

\description{
  This plots only the 25th to 75th percentile of expression intensities
  (Interquartile Range), sorted from smallest to largest IQR.  This
  modification is more readable than a regular boxplot for large sample
  sizes.  An optional batch variable may be specified, so that the sorting
  is done within each batch.  
}

\usage{
sortedIqrPlot(lumiobj, batchvar = rep(1, ncol(lumiobj)), ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{lumiobj}{
    An object of class LumiBatch, AffyBatch, ExpressionSet, or matrix.  This
    should contain raw, unnormalized, expression intensities.
}
\item{batchvar}{
  Optional integer batch variable.  If specified, samples will be sorted
  within batches only.  Default is to assume a single batch for all data.
}
\item{dolog2}{
  If TRUE, data will be log2-transformed before plotting.  Default is FALSE.
}
\item{\dots}{
  Optional arguments passed on to the plot() function.
}
}
\details{
  Function produces a plot only.
}
\value{
  none
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Levi Waldron
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library(ffpe)
data(lumibatch.GSE17565)
sortedIqrPlot(lumibatch.GSE17565,main="GSE17565")

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

