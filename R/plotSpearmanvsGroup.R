plotSpearmanvsGroup <-
function(lumibatch,
                                logtransform=TRUE,
                                goby=3,
                                xaxis=c("index","notindex"),
                                QCmeasure=c("IQR","ndetectedprobes","MAplot.var"),
                                cor.to=c("similar","pseudochip"),
                                pseudochip.samples=1:ncol(lumibatch),
                                detectionTh=0.01,
                                rankcutoff=NULL,
                                mincor=0,
                                maxcor=0.6,
                                below.smoothed.threshold=1.5,
                                lowess.f=1/3,
                                labelnote=NULL,
                                pch=1,lw=4,
                                linecol="red",
                                make.legend=TRUE,
                                main.title=NA,
                                ...){
  ##-----------------------------------------------------------------------------------------------
  ##-----------------------------------------------------------------------------------------------
  ## This function sorts samples by either interquartile range or number of probes called as detected.
  ## It then computes Spearman correlation within groups of samples which are similar according to this
  ## measure.  Correlations are plotted against either this measure or by rank, with a Loess curve.
  ## The resulting plot can be helpful in deciding how many samples to discard on the basis of IQR or
  ## number of detected probes.
  ##-----------------------------------------------------------------------------------------------
  ## VARIABLE DEFINITIONS
  ##-----------------------------------------------------------------------------------------------
  ##lumibatch: a lumibatch object.
  ##logtransform: if TRUE, values will be log2-transformed before calculating IQR.  If using only ndetectedprobes as a QC measure, this is irrelevant
  ##goby: the number of chips above and below the center chip to include when calculating rank correlations.  Ignored if cor.to="pseudochip".
  ##xaxis: Plot the actual QC measure (IQR or number of detected probes), or relative rank of the chips (1 being lowest)
  ##QCmeasure: IQR - use interquartile range for ranking arrays
  ##           ndetectedprobes - use number of detected probes for ranking arrays
  ##           MAplot.var - use variance of M in MA plot
  ##           c("IQR","ndetectedprobes") - use both
  ##           numeric vector with length equal to the number of samples - use this for ranking arrays
  ##cor.to: "similar" for correlation within a sliding window, "pseudochip" for correlation to the median pseudochip.
  ##labelnote: optional label for the plots, used only if QCmeasure is a numeric vector
  ##detectionTh: nominal detection p-value to consider a probe as detected or not (0.01 by default)
  ##pch: plotting character to be used for points (see ?par)
  ##lw: line width for Loess curve
  ## linecol: line color for Loess curve
  ## ...: other arguments passed on to plot()
  ##-----------------------------------------------------------------------------------------------
  ##-----------------------------------------------------------------------------------------------
  ##define the plotting function
  ##-----------------------------------------------------------------------------------------------
  ##-----------------------------------------------------------------------------------------------
  makeplots <- function(expr.dat,goby,xaxis,QC.measure,rankcutoff,mincor,below.smoothed.threshold,labelnote="QC",pch,lw,linecol,cor.to,pseudochip.samples,make.legend,...){
    if(any(is.na(QC.measure)|is.nan(QC.measure))){
      navals <- which(is.na(QC.measure)|is.nan(QC.measure))
      warning(paste(length(navals),"NA and NaN values of QC.measure were removed, along with the associated samples"))
      expr.dat <- expr.dat[,-navals]
      QC.measure <- QC.measure[-navals]
    }
    expr.sort <- expr.dat[,order(QC.measure)]
    pseudochip.samples <- pseudochip.samples[order(QC.measure)]
    if(cor.to[1]=="similar"){
      if(is.na(main.title)|is.null(main.title)){
        thismain <- paste("Rank correlation of groups of ",2*goby+1," arrays\n grouped by ",labelnote,", lowest to highest",sep="")
      }else{
        thismain <- main.title
      }
      print("Calculating correlation within sliding window...")
      cor.vector <- rep(NA,length(QC.measure))
      for (i in (1+goby):(length(QC.measure)-goby)){
        thiscor <- cor(expr.sort[,(i-goby):min((i+goby),ncol(expr.sort))],method="spearman",use="pairwise.complete.obs")
        cor.vector[i] <- median(thiscor[-(1+goby),1+goby],na.rm=TRUE)
      }
      ##For the first few samples, calculate the median correlation to the first goby samples
      thiscor <- cor(expr.sort[,1:(1+goby)],method="spearman",use="pairwise.complete.obs")
      for (i in 1:(1+goby)){
        cor.vector[i] <- median(thiscor[-i,i],na.rm=TRUE)
      }
      names(cor.vector) <- colnames(expr.sort)
      for (i in 1:length(cor.vector)){
        if (is.na(cor.vector[i])){  ## replace NA values
          if(i < (length(cor.vector)-goby)){  #for the low QC values, NA values should be replaced with 0
            cor.vector[i] <- 0
          }else{
            ##otherwise, replace the last few NA values with the mean correlation of 5 arrays at the top of the QC range.
            cor.vector[i] <- mean(cor.vector[(length(cor.vector)-5):length(cor.vector)],na.rm=TRUE) 
          }
        }
      }
    }else if(cor.to[1]=="pseudochip"){
      if(is.na(main.title)|is.null(main.title)){
        thismain <- paste("Rank correlation to median pseudochip \n grouped by ",labelnote,", lowest to highest",sep="")
      }else{
        thismain <- main.title
      }
      print("Calculating Spearman correlation to median pseudochip using pairwise complete observations.")
      print(paste("Using samples:",paste(colnames(expr.sort)[pseudochip.samples],collapse=", ")))
      pseudochip <- apply(expr.sort[,pseudochip.samples],1,median,na.rm=TRUE)
      cor.vector <- cor(expr.sort,pseudochip,method="spearman",use="pairwise.complete.obs")[,1]
    }
    QC.measure.sort <- QC.measure[order(QC.measure)]
    if(xaxis[1]=="index"){
      cormat.2col <- data.frame(i=1:length(cor.vector),
                                spearman=cor.vector,
                                row.names=names(cor.vector))
    }else{
      cormat.2col <- data.frame(i=QC.measure.sort,
                                spearman=cor.vector,
                                row.names=names(cor.vector))
    }
    cormat.2col <- cormat.2col[!is.na(cormat.2col$spearman),]
    cor.vector <- cor.vector[match(rownames(cormat.2col),names(cor.vector))]
    rejectQC <- rep(FALSE,length(cor.vector));names(rejectQC) <- names(cor.vector)
    if(!is.na(lowess.f)&!is.null(lowess.f)){
      movingavg.n <- 2*goby+1
      library(TTR)
myavg <- rev(SMA(rev(cormat.2col$spearman),n=movingavg.n))
      if(class(myavg)=="try-error"){
        warning(paste("error in loess fit for",labelnote))
        if(is.null(rankcutoff)) rankcutoff <- 1
      }else{
        was.na <- as.integer(attributes(na.omit(as.numeric(myavg)))$na.action)
        ##Use original values for those that couldn't be smoothed by SMA:
        myavg[was.na] <- cormat.2col[was.na,"spearman"]
        cormat.2col$movingaverage <- myavg
        cormat.loess <- loess(movingaverage~i,data=cormat.2col,span=lowess.f,)
        cormat.2col$interpolate.i <- seq(min(cormat.2col$i),max(cormat.2col$i),length.out=nrow(cormat.2col))
        cormat.2col$smoothed <- predict(cormat.loess,data.frame(i=cormat.2col$interpolate.i))
        library(sfsmisc)
        cormat.2col$ddy <- D1D2(cormat.2col$interpolate.i,cormat.2col$smoothed,deriv=2,spar.offset=0.6)$D2
        cormat.2col$ddy[1] <- NA
        ##only consider samples with sufficiently low correlation as the cutoff:
        cormat.lowspearman <- cormat.2col[cormat.2col$spearman<maxcor,]
        if(is.null(rankcutoff)) rankcutoff <- cormat.lowspearman$interpolate.i[which.min(cormat.lowspearman$ddy)]
        diff.from.smoothed <- cormat.2col$spearman-cormat.2col$smoothed
        reject.below.smoothed <- diff.from.smoothed < (-below.smoothed.threshold*IQR(diff.from.smoothed))
      }
      if(xaxis[1]!="index"){
        rankcutoff <- sum(cormat.2col$i < rankcutoff)
      }
      rejectQC[1:rankcutoff] <- TRUE
      if(exists("reject.below.smoothed")) rejectQC[reject.below.smoothed] <- TRUE
    } #end if(!is.na(lowess.f)&!is.null(lowess.f))
    rejectQC[cor.vector<mincor] <- TRUE
    if(all.equal(names(rejectQC),rownames(cormat.2col))){
      cormat.2col$rejectQC <- rejectQC
    }else{
      warning("Something is wrong with the sample names.")
    }
    mycol <- ifelse(rejectQC,"red","black")
    ##correct label note
    if(xaxis[1]=="index"){
      xlab <- paste(labelnote,"rank, 1 is smallest")
    }else{
      xlab <- labelnote
    }
    ##make the plot:
    plot(spearman~i,data=cormat.2col,
         pch=pch,
         col=mycol,
         xlab=xlab,
         ylab="Spearman correlation",
         main=thismain,
         ...)
    if(mincor>0) abline(h=mincor,col="red")
    abline(v=cormat.2col[rankcutoff,"i"],col="red")
    reject.summary <- summary(cormat.2col$rejectQC)
    if(make.legend){
      legend("bottomright",pch=1,lty=-1,legend=c(paste(reject.summary[3],"rejected"),paste(reject.summary[2],"not rejected")),bty='n',col=c("red","black"))
    }
    if("smoothed" %in% colnames(cormat.2col)){
      lines(smoothed~interpolate.i,data=cormat.2col,lw=lw,col="red")
    }
    return(cormat.2col)
  }
  ## get ready then call the plotting function
  if(logtransform){
    expr.dat <- log2(exprs(lumibatch))
  }else{
    expr.dat <- exprs(lumibatch)
  }
  output <- list()
  if(identical(class(QCmeasure),"numeric") | identical(class(QCmeasure),"integer")){
    for (thisxaxis in xaxis){
      if(is.na(labelnote)|is.null(labelnote)) labelnote <- "custom QC"
      thismethod <- paste(thisxaxis,QCmeasure,sep="_")
      output[[thismethod]] <- makeplots(expr.dat,goby,thisxaxis,QC.measure=QCmeasure,cor.to=cor.to,pseudochip.samples=pseudochip.samples,rankcutoff,mincor,below.smoothed.threshold=below.smoothed.threshold,labelnote,lw=lw,pch=pch,linecol=linecol,make.legend=make.legend,...)
    }
  }else if("IQR"%in% QCmeasure){
    Raw.IQR <- apply(expr.dat,2,IQR,na.rm=TRUE)
    Raw.IQR[is.na(Raw.IQR)] <- 0
    for (thisxaxis in xaxis){
      thismethod <- paste(thisxaxis,QCmeasure,sep="_")
      output[[thismethod]] <- makeplots(expr.dat,goby,thisxaxis,QC.measure=Raw.IQR,cor.to=cor.to,pseudochip.samples=pseudochip.samples,rankcutoff,mincor,below.smoothed.threshold=below.smoothed.threshold,labelnote="IQR",lw=lw,pch=pch,linecol=linecol,make.legend=make.legend,...)
    }
  }else if("ndetectedprobes"%in% QCmeasure){
    Raw.fractiondetected <- apply(detection(lumibatch),2,function(x) sum(x<detectionTh)/length(x))
    for (thisxaxis in xaxis){
      thismethod <- paste(thisxaxis,QCmeasure,sep="_")
      output[[thismethod]] <- makeplots(expr.dat,goby,thisxaxis,QC.measure=Raw.fractiondetected,cor.to=cor.to,pseudochip.samples=pseudochip.samples,rankcutoff,mincor,below.smoothed.threshold=below.smoothed.threshold,labelnote="fraction detected",lw=lw,pch=pch,linecol=linecol,make.legend=make.legend,...)
    }
  }else if("MAplot.var"%in% QCmeasure){
    ma.results <- ma.stats(expr.dat)
    M <- na.omit(ma.results[["M"]])
    M.var <- 1/apply(M,2,function(x){  ##variance of difference from Lowess curve
      output <- try(var(x-lowess(x)$y))
      if(class(output)=="try-error") output <- 1000
      return(output)
    })
##    M.var <- 1/apply(M,2,function(x) sum(abs(x))/length(x))
##    M.var <- 1/apply(M,2,var)
    M.var[is.na(M.var)] <- 0
    for (thisxaxis in xaxis){
      thismethod <- paste(thisxaxis,QCmeasure,sep="_")
      output[[thismethod]] <- makeplots(expr.dat,goby,thisxaxis,QC.measure=M.var,cor.to=cor.to,pseudochip.samples=pseudochip.samples,rankcutoff,mincor,below.smoothed.threshold=below.smoothed.threshold,labelnote="IQR",lw=lw,pch=pch,linecol=linecol,make.legend=make.legend,...)
    }
  }
  ##Sort the output to the original sample order
  output <- lapply(output,function(x) x[na.omit(match(sampleNames(lumibatch),rownames(x))),])
  if(length(output)==1) output <- output[[1]]
  return(output)
}

