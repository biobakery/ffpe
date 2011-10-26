sortedIqrPlot <-
function(data.obj,batchvar=rep(1,ncol(data.obj)),dolog2=FALSE,...){
  if(!length(batchvar)==ncol(data.obj)) stop("length(batchvar) must be equal to ncol(data.obj)")
  if(class(data.obj)=="LumiBatch")
    {
      library(lumi)
      data.obj <- exprs(data.obj)
    }else if(class(data.obj)=="AffyBatch" | class(data.obj)=="ExpressionSet"){
      library(affy)
      data.obj <- exprs(data.obj)
    }else if(!class(data.obj)=="matrix"){
      stop("data.obj should be of class LumiBatch, AffyBatch, ExpressionSet, or matrix.")
    }
  if(dolog2)
    {
      if(min(data.obj) <= 0)
        {
          data.obj <- data.obj - min(data.obj) + 1
        }
      data.obj <- log2(data.obj)
    }
  persample.iqr <- apply(data.obj,2,quantile,probs=c(0.25,0.75),na.rm=TRUE)
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

