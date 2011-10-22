sortedIqrPlot <-
function(lumiobj,batchvar=rep(1,ncol(lumiobj)),dolog2=FALSE,...){
  if(!length(batchvar)==ncol(lumiobj)) stop("length(batchvar) must be equal to ncol(lumiobj)")
  if(class(lumiobj)=="LumiBatch")
    {
      library(lumi)
      lumiobj <- exprs(lumiobj)
    }else if(class(lumiobj)=="AffyBatch" | class(lumiobj)=="ExpressionSet"){
      library(affy)
      lumiobj <- exprs(lumiobj)
    }else if(!class(lumiobj)=="matrix"){
      stop("lumiobj should be of class LumiBatch, AffyBatch, ExpressionSet, or matrix.")
    }
  if(dolog2)
    {
      if(min(lumiobj) <= 0)
        {
          lumiobj <- lumiobj - min(lumiobj) + 1
        }
      lumiobj <- log2(lumiobj)
    }
  persample.iqr <- apply(lumiobj,2,quantile,probs=c(0.25,0.75),na.rm=TRUE)
  persample.iqr <- persample.iqr[,order(batchvar,apply(persample.iqr,2,diff))]
  plot(1:ncol(persample.iqr),
       seq(min(persample.iqr),max(persample.iqr),length.out=ncol(persample.iqr)),
       type='n',
       ...)
  for (i in 1:ncol(persample.iqr))
    {
      segments(x0=i,
               y0=persample.iqr[1,i],
               y1=persample.iqr[2,i])
    }
}

