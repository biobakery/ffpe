\name{CATplot}
\alias{CATplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CATplot(vec1, vec2, maxrank, make.plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec1}{
%%     ~~Describe \code{vec1} here~~
}
  \item{vec2}{
%%     ~~Describe \code{vec2} here~~
}
  \item{maxrank}{
%%     ~~Describe \code{maxrank} here~~
}
  \item{make.plot}{
%%     ~~Describe \code{make.plot} here~~
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
function(vec1,vec2,maxrank,make.plot=TRUE,...){ 
  ##----------------------------------------------------------------------
  ##----------------------------------------------------------------------
  ##Fuction for Concordance at the Top plots
  ##----------------------------------------------------------------------
  ## vec1 and vec2: ranked lists to use for CAT plot
  ## maxrank: only go up to maxrank
  ## make.plot: if TRUE, the plot will be made
  ## ...: arguments passed on to plot
  ##----------------------------------------------------------------------  
  x1 <- sort(vec1)
  x2 <- sort(vec2)
  if(is.na(maxrank)|is.null(maxrank)|maxrank>min(length(x1),length(x2))) maxrank <- min(length(vec1),length(vec2))
  output <- data.frame(rank=1:maxrank,
                       fraction=(1:maxrank)/max(length(x1),length(x2)),
                       concordance=NA)
  for (i in 1:nrow(output)){
    output[i,"concordance"] <- length(intersect(names(x1)[1:i],names(x2)[1:i]))/i
  }
  if(make.plot) plot(concordance~rank,data=output,type='l',...)
  return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
