\name{sortedIqrPlot}
\alias{sortedIqrPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sortedIqrPlot(lumiobj, batchvar = rep(1, ncol(lumiobj)), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lumiobj}{
%%     ~~Describe \code{lumiobj} here~~
}
  \item{batchvar}{
%%     ~~Describe \code{batchvar} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
function(lumiobj,batchvar=rep(1,ncol(lumiobj)),...){
  ##----------------------------------------------------------------------
  ##----------------------------------------------------------------------
  ##minimal box plot of samples, which get ordered optionally by batch then IQR.
  ##----------------------------------------------------------------------
  ## lumiobj: lumi object or expression matrix, should be log2-transformed already
  ## batchvar: optional integer batch variable.  If specified, samples will be sorted within batches only.
  ##----------------------------------------------------------------------
  if(!length(batchvar)==ncol(lumiobj)) stop("length(batchvar) must be equal to ncol(lumiobj)")
  if(class(lumiobj)=="LumiBatch" | class(lumiobj)=="expressionSet") lumiobj <- exprs(lumiobj)
  persample.iqr <- apply(lumiobj,2,quantile,probs=c(0.25,0.75),na.rm=TRUE)
  persample.iqr <- persample.iqr[,order(batchvar,apply(persample.iqr,2,diff))]
  plot(1:ncol(persample.iqr),seq(min(persample.iqr),max(persample.iqr),length.out=ncol(persample.iqr)),type='n',...)
  for (i in 1:ncol(persample.iqr)) segments(x0=i,y0=persample.iqr[1,i],y1=persample.iqr[2,i])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
