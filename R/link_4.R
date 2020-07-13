## arguments:
## k.mcmc   : K output as three dimensional array from BDgraph
## x   : index vector of source nodes to link
## y   : index vector of target nodes to link
## delta    : tolerence thresholds for "practical significance"
##            just in case many links are identified?
## alpha : probability for significance
## partial: whether to compute partial correlation
## names : the names of the variables (default extracted from the obj.)
## expand : whether to expand x and y or treat them as is
## flag   : flag to denote how to calculate the correlation graph
##          1 use inverse of precision matrix; 2 use the latent variable
##          3 use rank correlation of the predictive distribution
## value    : posterior estimates of the partial correlation
##            coefficients for the links
link <- function(k.mcmc,x,y=NULL,delta=0.05,alpha=0.05,partial=TRUE,
                 names=NULL,expand=FALSE,flag=1)
{
  if(any(!(x %in% 1:dim(k.mcmc)[2]))){
    stop("x misspecified.\n ")
  }
  if(any(!(y %in% 1:dim(k.mcmc)[2]))){
    stop("y misspecified.\n ")
  }
  if(is.null(names)){
    names <- as.character(seq(1,nrow(k.mcmc[1,,])))
  }
  ## storage
  rho <- array(NA,dim(k.mcmc))
  if(partial){
    for(i in 1:dim(k.mcmc)[1]){
      ## loop through the MC iterations
      ## standardize the precision matrix
      #theta <- solve(cov2cor(solve(k.mcmc[i,,])))
      theta <- k.mcmc[i,,] ## no standardization?
      ## maintain the sparsity of the standardized precision matrix
      #theta[abs(k.mcmc[i,,])<1e-6] <- 0
      ## compute the partial correlation coefficients
      rho[i,,] <- -1*cov2cor(theta)
    }
  }
  else{
    for(i in 1:dim(k.mcmc)[1]){
      cat("posterior sampe ",i,"\n")
      ## loop through the MC iterations
      ## compute the marginal correlation coefficients
      if(flag==1){
        ## inverse of precision matrix
        rho[i,,] <- cov2cor(solve(k.mcmc[i,,]))
      }
      # else if (flag==2){
      #   ## correlation of latent variables
      #   rho[i,,] <- cor(mcmc$all_Z[i,,])
      # }
      # else if (flag==3){
      #   ## rank correlation of predictive distribution
      #   rho[i,,] <- cor(mcmc$all_Y[i,,],method="kendall")
      # }
    }
  }
  ## return
  if(expand){
    ## y ignored
    result <- expand.grid(source=x,target=x)
    result <- subset(result,source > target)
  }
  else if(is.null(y)){
    result <- expand.grid(source=x,target=x)
    result <- subset(result,source > target)
  }
  else{
    result <- data.frame(source=x,target=y)
  }
  
  Ests <- vector("list",nrow(result))
  ## posterior correlation coefficients between source and target
  for(j in 1:nrow(result)){
    ## loop through the pairs
    l.source <- result$source[j]
    l.target <- result$target[j]
    ## posterior estimates and uncertainty
    l.q <- quantile(rho[,l.source,l.target],probs=c(0.5,alpha/2,1-alpha/2))
    names(l.q) <- NULL
    l.p <- signif4(rho[,l.source,l.target],delta = delta)
    names(l.p) <- c("negative","zero","positive","ratio","probSig")
    ## determine the stength of association
    Ests[[j]] <- c(est=l.q[1],lwr=l.q[2],upr=l.q[3],
                   sgnf=signif1(l.q[2],l.q[3],l.q[1],thresh=0.0001),
                   l.p)
  }
  #browser()
  result$source <- names[result$source]
  result$target <- names[result$target]
  cbind(result,do.call(rbind,Ests))
}

signif1 <- function(lwr,upr,est,thresh=0.05){
  c1 <- (round(lwr,2)*round(upr,2))>=0
  c2 <- abs(est)>thresh
  c1 & c2
}

signif4 <- function(x,delta=0.05){
  ## x: coefficients samples
  ## value: probability of practically significant
  ## defined as the likelihood ratio of the dominant portion
  ## over the rest
  p1 <- mean(x>delta) ## exceeds positive indifference
  p2 <- mean(x<(-1*delta)) ## below negative indifference
  p0 <- mean(abs(x)<delta) ## indifference (not practically important)
  
  if(p1>=p2){ ## most mass in negative zone
    ll <- p1/(p0+p2)
    pp <- p1
  }else{
    ## most mass in positive zone
    ll <- p2/(p0+p1)
    pp <- p2
  }
  c(p2,p0,p1,ll,pp)
}

test_signf <- function(){
  rm(list=ls())
  source("link_4.R")
  signif4(rnorm(100))
  signif4(rnorm(100,sd=0.1))
  signif4(rnorm(100,sd=10))
  signif4(rnorm(100,mean=0.1))
  signif4(rnorm(100,mean=0.1,sd=0.1))
}

signif2 <- function(lwr,upr,thresh=0){
  c1 <- (lwr * upr) > 0 ## significant
  lwr1 <- ifelse(abs(lwr)<thresh,0,lwr)
  upr1 <- ifelse(abs(upr)<thresh,0,upr)
  c2 <- (lwr1 * upr1) >=0 ## approximately significant
  c1 | c2
}

