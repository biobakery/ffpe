sortedIqrPlot <-
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

