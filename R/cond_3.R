## cond_3.R
## posterior predictive check based conditional probability
## x: posterior sample object from BDgraph
## formula: define the recruitment response 
##         and environmental variable to run posterior predictive check
## data: observed data to evaluate the formula
## nquantile: number of quantiles to evaluate the ECDF of response
## breaks : the break point to discretize the environmental variable
##          default (below and above arithmetic mean)
## label  : the label to define each discretized environmental variable
##          default (below and above)
## value  : an 3 dimensional array
##        first dimension is the posterior samples
##        second dimension is the discritized response variables
##        third  dimension is the discretized environmental variable
##        value is the cumulative probability of response below certain value
##        value is missing if no environmental variable satisfy the specified condition.
##         average 
cond <- function(x,formula,data,breaks=NULL,labels=NULL,nquantile=49)
{
  v <- all.vars(formula)
  stopifnot(all(v %in% names(data)))
  
  ## browser()
  ## the response variable grid 
  l.y <- data[,v[1]]
  l.y.grid <- quantile(
    l.y,probs = seq(0,1,length.out = nquantile),
    na.rm=TRUE)
  
  ## discretize the environmental variable
  l.x <- data[,v[2]]
  if(is.null(breaks)){
    breaks <- c(-Inf,mean(l.x,na.rm=T),Inf)
    labels <- c("Below Average","Above Average")
  }
  
  ## posterior predictive distributions
  l.pred.ii <- match(v,row.names(x$last_graph))
  l.pred.y <- x$all_Y[,,l.pred.ii[1],drop=TRUE] ## response
  l.pred.x <- x$all_Y[,,l.pred.ii[2],drop=TRUE] ## environmental variable
  
  l.result <- array(
    NA,c(nrow(l.pred.y),
         length(l.y.grid),
         length(labels))
    )
  dimnames(l.result) <- list(
    NULL,
    paste(v[1],names(l.y.grid),sep=" "),
    paste(v[2],labels,sep=" ")
  )
  for(i in 1:nrow(l.pred.y)){ ## for each posterior sample
    ## discretize the environmental variable
    l.x.cut <- cut(l.pred.x[i,],breaks=breaks,labels=labels)
    ## divide the time series according to environmental condition
    l.y.lst <- split(l.pred.y[i,],f = l.x.cut)
    ## compute the ECDF for each group of time periods
    l.r.lst <- vector("list",length(l.y.lst))
    for(j in 1:length(l.y.lst)){
      ## ECDF of reponse among years satisfying environmental condition
      l.Fn <- ecdf(l.y.lst[[j]]) 
      ## store ECDF
      l.r.lst[[j]] <- l.Fn(l.y.grid)
    }
    ## store as matrix
    l.r <- matrix(NA,length(l.y.grid),length(labels))
    l.r[,match(names(l.y.lst),labels)] <- do.call(cbind,l.r.lst)
    ## store as array
    l.result[i,,] <- l.r
  }
  result <- list(value=l.result,dummy_=l.y.grid,formula=formula,
                 breaks=breaks,labels=labels)
  names(result)[2] <- v[1]
  result
}

## helper function to visualize the results
## x: output from conditional probability
## col: vector of colors to visualize each environmental category
## lwd: vector of line widths to visualize each category
## fill: vector of colors to fill each environmental category
## alpha: controls the posterior credible levels
## title: default the environmental variable
## title.position: legend poisition for title
## text.position: position for text 
## alpha.f: fill transparency 0.6
plot.cond <- function(
  x,col=c("white","white"),fill=c("salmon","steelblue1"),
  lty=c(1,1),lwd=4,
  alpha=0.1,title=NULL,
  title.position=NULL,text.position=NULL,
  alpha.f=0.6,ylim=c(-0.5,0.5),...){
  #browser()
  ecdf.diff <- x$value[,,2]-x$value[,,1]
  ecdf.median <- apply(ecdf.diff,2,median)
  ecdf.lwr <- apply(ecdf.diff,2,quantile,probs=alpha/2)
  ecdf.upr <- apply(ecdf.diff,2,quantile,probs=1-alpha/2)
  xlim__ <- range(x[[2]])
  ylim__ <- ylim
  plot(x[[2]],ecdf.median,xlim=xlim__,ylim=ylim__,
          type="l",col=col[1],lwd=lwd,lty=lty[1],...)
  #l.fill <- cbind(fill,rep(1,ncol(ecdf.lwr)))[,1]
  #for(i in 1:ncol(ecdf.lwr)){
  l.poly <- cbind(c(x[[2]],rev(x[[2]])),c(ecdf.lwr,rev(ecdf.upr)))
  polygon(l.poly,col=adjustcolor(fill[1],alpha.f=alpha.f),
            border="black",lty=2,lwd=0.5)
  
  lines(x[[2]],ecdf.median,xlim=xlim__,ylim=ylim__,
          type="l",col=col[1],lwd=lwd,lty=lty[1])
  #matlines(x[[2]],ecdf.median,xlim=xlim__,ylim=ylim__,
  #         type="l",col=l.fill,lwd=1,lty=lty)
  
  if(is.null(title)){
    title <- names(x)[2]
  }
  if(is.null(title.position)){
    title.position <- c(x[[2]][2],0.9)
  }
  text(title.position[1],title.position[2],labels = title,adj=0)
  
  delta.ecdf <- apply(x$value,1,function(mat){
    mean(mat[,2]-mat[,1])
  })
  delta <- quantile(delta.ecdf,probs=c(alpha/2,1-alpha/2))
  text__ <- paste(round(delta,2),collapse=",")
  if(is.null(text.position)){
    text.position <- c(x[[2]][length(x[[2]])],0.3)
  }
  text(text.position[1],text.position[2],
       labels=bquote(Delta(p)~(.(text__))),adj=1)
  grid()
}
test <- function(){
  rm(list=ls())
  load(file="pnm3m.rda")
  load(file="pnm_result_3m.rda")
  source("cond_2.R")
  #out <- cond(x=sample.mcmc,formula=age0~NAO_Jun,data=pnm2)
  out <- cond(x=sample.mcmc,formula=age0~DO,data=pnm2)
  debug(plot.cond)
  plot.cond(out,xaxt="n",xlab="",ylab="")
  # ecdf.mean <- apply(out,c(2,3),mean)
  # plot(ts(ecdf.mean),plot.type="single",
  #      lwd=4,lty=c(2,1),col=c(grey(0),grey(0.5)))
  # delta.ecdf <- apply(out,1,function(mat){
  #   delta <- mat[,1]-mat[,2]
  #   mean(delta)
  # })
  # hist(delta.ecdf)
  # delta.ecdf <- apply(out,1,function(mat){
  #   mat[,1]-mat[,2]
  # })
  # delta.ecdf.lmt <- t(apply(delta.ecdf,1,quantile,prob=c(0.05,0.5,0.95)))
  # plot(ts(delta.ecdf.lmt),plot.type="single",lty=c(2,1,2),lwd=c(1,2,1))
  # abline(h=0,lwd=2,col="red")
  # 
  # ecdf.median <- apply(out,c(2,3),median)
  # ecdf.lwr <- apply(out,c(2,3),quantile,prob=0.25)
  # ecdf.upr <- apply(out,c(2,3),quantile,prob=0.75)
  # plot(ts(cbind(ecdf.median,ecdf.lwr,ecdf.upr)),
  #      plot.type="single",
  #      lty=c(1,1,2,2,2,2),lwd=c(4,4,1,1,1,1),
  #      col=rep(c(grey(0),grey(0.5)),each=3))
  
}