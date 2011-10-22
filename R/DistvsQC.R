DistvsQC <-
function(mat,Raw.IQR,...){
  ##----------------------------------------------------------------------
  ##----------------------------------------------------------------------
  ## This function calculates pairwise sample correlation as a function of the
  ## difference, mean, minimum, and geometric mean of some continuous QC measure of the chips
  ##----------------------------------------------------------------------
  ## mat: expression matrix
  ## Raw.IQR: IQR of raw expression data, or another continuous measure available for each chip
  ## ... are arguments passed on to cor()
  ##----------------------------------------------------------------------
  distmat <- 1-cor(mat)
  IQRmean=outer(X=Raw.IQR,Y=Raw.IQR,FUN="+")/2
  IQRdiff=outer(X=Raw.IQR,Y=Raw.IQR,FUN="-")
  IQRgeomean=sqrt(outer(X=Raw.IQR,Y=Raw.IQR,FUN="*"))
  IQRmin=outer(X=Raw.IQR,Y=Raw.IQR,FUN="pmin")
  all.equal(names(Raw.IQR),colnames(distmat))
  dist.IQR.df <- data.frame(dist=distmat[upper.tri(distmat)],
                            IQRmean=IQRmean[upper.tri(IQRmean)],
                            IQRgeomean=IQRgeomean[upper.tri(IQRgeomean)],
                            IQRmin=IQRmin[upper.tri(IQRmin)],
                            IQRdiff=IQRdiff[upper.tri(IQRdiff)])
  convertCutLabels <- function(myfac){
    tmp <- levels(myfac)
    tmp <- sub("(","",tmp,fixed=TRUE)
    tmp <- sub("]","",tmp,fixed=TRUE)
    tmp <- strsplit(tmp,",")
    tmp <- sapply(tmp,function(x) mean(as.numeric(x)))
    tmp <- as.character(tmp)
    levels(myfac) <- tmp
    return(myfac)
  }
  dist.IQR.df$IQRmindiscrete <- cut(dist.IQR.df$IQRmin,breaks=seq(from=0,to=max(dist.IQR.df$IQRmin),by=0.1))
  dist.IQR.df$IQRmindiscrete <- convertCutLabels(dist.IQR.df$IQRmindiscrete)
  dist.IQR.df$IQRgeomeandiscrete <- cut(dist.IQR.df$IQRgeomean,breaks=seq(from=0,to=max(dist.IQR.df$IQRgeomean),by=0.1))
  dist.IQR.df$IQRgeomeandiscrete <- convertCutLabels(dist.IQR.df$IQRgeomeandiscrete)
  return(dist.IQR.df)
}

