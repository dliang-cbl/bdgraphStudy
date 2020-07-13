## DIC computation for Gaussian Coupla model
dic_calc <- function(x,z,index=NULL,marginal=FALSE)
{
  #n <- nrow(z)
  #S <- var(z)*(n-1)
  meanK <- apply(x$all_K,c(2,3),mean)
  Dhat <- -2*log.lik(meanK,z,index=index,marginal=marginal)
  Dbar.mat <- apply(x$all_K,1,function(k){
    -2*log.lik(k,z,index=index,marginal=marginal)
  })
  Dbar <- mean(Dbar.mat)
  c(Dbar=Dbar,pD=Dbar-Dhat,DIC=2*Dbar-Dhat)
}
coef_ggm <- function(x,z,index,xvar=NULL){
  ##x: fitted ggm object
  ##z: observed raw data
  ## index: response variable (must be univariate)
  ## xvar: predictors to group together (default is everything else)
  ## value: fitted mean and uncertainty given the dat
  
  ## residuals (local storage)
  r <- scale(z,scale = FALSE)
  ## marginal mean
  mu.y <- apply(z[,index,drop=FALSE],2,mean)
  
  if(is.null(xvar)){
    xvar <- -1* index
  }
  store <- vector("list",dim(x$all_K)[1])
  for(i in 1:dim(x$all_K)[1]){
    K <- x$all_K[i,,]
    
    ## standardize variates
    rho <- cov2cor(solve(K)) ## correlation matrix
    theta <- solve(rho)  ## precision matrix after standardization
    
    ## regression mean
    beta <- -1*theta[index,xvar]/theta[index,index]
    coef <- beta %*% t(r[,xvar,drop=FALSE])

    store[[i]] <- c(coef) 
  }
  
  tmp <- abind(store,along=-1)
  t(apply(tmp,2,function(x){
    c(mean(x),sd(x))
  }))
}

## simulation based estimates of confiditional mean and variance
fitted_ggm <- function(x,z,index=NULL,type="pred"){
  ## residuals (local storage)
  r <- scale(z,scale = FALSE)
  ## marginal mean
  mu.y <- apply(z[,index,drop=FALSE],2,mean)
  
  store <- vector("list",dim(x$all_K)[1])
  for(i in 1:dim(x$all_K)[1]){
    #K <- apply(x$all_K,c(2,3),mean)
    K <- x$all_K[i,,]
    
    ## standardize variates
    rho <- cov2cor(solve(K)) ## correlation matrix
    theta <- solve(rho)  ## precision matrix after standardization
    
    ## regression mean
    a0 <- solve(rho[-index,-index,drop=FALSE],t(r[,-index,drop=FALSE]))
    a1 <- rho[index,-index,drop=FALSE] %*% a0
    mu.cond.y <- mu.y + t(a1)
    ## regression variance
    sigma.cond.y <- rho[index,index,drop=FALSE] - rho[index,-index,drop=FALSE] %*% solve(
      rho[-index,-index,drop=FALSE],rho[-index,index,drop=FALSE]
    )
    ## this won't work for multiple variables
    if(type=="conf"){
      store[[i]] <- c(mu.cond.y) #uncertainty of the mean only
    }else{
      store[[i]] <- rnorm(n=nrow(mu.cond.y),mean=mu.cond.y,sd= sqrt(sigma.cond.y))  
    }
  }
  #browser()
  tmp <- abind(store,along=-1)
  t(apply(tmp,2,function(x){
    c(mean(x),sd(x))
  }))
}
log.lik <- function(K,z,index=NULL,marginal=FALSE,loo=FALSE){
  ## K: precision matrix
  ## z: observed data
  ## n: sample size
  ## index: the subset of variables to focus (default ALL)
  ## marginal: whether compute the marginal (TRUE) or conditional (FALSE) density
  ## loo: logical, whether remove a vector for each individual observations
  ##      or a single value
  ## value: log-likelihood.
  ## browser()
  myvar <- function(y,mu=NULL){
    ## calculate empirical covariance
    p <- ncol(y)
    if(is.null(mu)){
      mu <- apply(y,2,mean)
    }
    r <- sweep(y,2,mu)
    tmp <- apply(r,1,function(x) outer(x,x,FUN="*"))
    r2 <- apply(tmp,1,sum)
    dim(r2) <- c(p,p)
    r2
  }
  ## standardize variates
  rho <- cov2cor(solve(K)) ## correlation matrix
  theta <- solve(rho)  ## precision matrix after standardization
  ## residuals (local storage)
  r <- scale(z,scale = FALSE)
  if(is.null(index)){
    n <- nrow(z)
    S <- var(z)*(n-1)
    prod <- crossprod(theta,S)
    if(loo){
      llik <- log.lik.loo(theta,z)
    }else{
      llik <- 0.5*n*log(det(theta)) - 0.5*sum(diag(prod))
    }
  }
  else{
    theta0 <- matrix(NA,length(index),length(index)) ## storage for precision matrix
    ## has focus variables
    if(marginal){ ## marginal density
      n <- nrow(z[,index,drop=FALSE])
      S <- var(z[,index,drop=FALSE])*(n-1)
      theta0 <- solve(rho[index,index,drop=FALSE])
      prod <- crossprod(theta0,S)
      if(loo){
        llik <- log.lik.loo(theta0,z[,index,drop=FALSE])
      }else{
        llik <- 0.5*n*log(det(theta0)) - 0.5*sum(diag(prod))
      }
    }else{ ## conditional density of focus variable given the rest
      ## marginal mean
      mu.y <- apply(z[,index,drop=FALSE],2,mean)
      ## regression mean
      a0 <- solve(rho[-index,-index,drop=FALSE],t(r[,-index,drop=FALSE]))
      a1 <- rho[index,-index,drop=FALSE] %*% a0
      mu.cond.y <- mu.y + t(a1)
      ## regression variance
      sigma.cond.y <- rho[index,index,drop=FALSE] - rho[index,-index,drop=FALSE] %*% solve(
        rho[-index,-index,drop=FALSE],rho[-index,index,drop=FALSE]
      )
      ## regression precision
      theta0 <- solve(sigma.cond.y)
      ## SSE
      n <- nrow(z)
      a2 <- z[,index,drop=FALSE]-mu.cond.y
      a3 <- apply(a2,1,function(x) outer(x,x,FUN="*"))
      if(is.null(dim(a3))){
        S <- matrix(sum(a3),length(index),length(index))
      }else{
        S <- matrix(apply(a3,1,sum),length(index),length(index))
      }
      #S <- var(z[,index]-mu.cond.y)*(n-1)
      prod <- crossprod(theta0,S)
      if(loo){
        llik <- log.lik.loo(theta0,z[,index,drop=FALSE],mu=mu.cond.y)
      }else{
        llik <- 0.5*n*log(det(theta0)) - 0.5*sum(diag(prod))
      }
    }
  }
  llik
}

test <- function()
{
  rm(list=ls())
  load(file="pnm.rda")
  load(file="pnm_result.rda")
  copula.transfer <- function(data){
    FinvData <- apply(data,2,rank,ties.method="random")/(nrow(data)+1)
    qnorm(FinvData)
  }
  source("dic_3.R")
  set.seed(123)
  Ztilde <- copula.transfer(pnm2)
  log.lik.loo(theta=sample.mcmc$all_K[1,,],Ztilde=Ztilde)
  
  log.lik(sample.mcmc$all_K[1,,],z=Ztilde,index=1:20,marginal=TRUE)
  sum(log.lik(sample.mcmc$all_K[1,,],z=Ztilde,index=1:20,marginal=TRUE,loo=T))
  ## -167.4875
  
  log.lik(sample.mcmc$all_K[1,,],z=Ztilde,index=21:41,marginal=FALSE)
  sum(log.lik(sample.mcmc$all_K[1,,],z=Ztilde,index=21:41,marginal=FALSE,loo=T))
  ## -139.9296
  
  log.lik(sample.mcmc$all_K[1,,],z=Ztilde,index=1,marginal=TRUE)
  sum(log.lik(sample.mcmc$all_K[1,,],z=Ztilde,index=1,marginal=TRUE,loo=T))
  ## -11.73
  
  log.lik(sample.mcmc$all_K[1,,],z=Ztilde,index=2:41,marginal=FALSE)
  sum(log.lik(sample.mcmc$all_K[1,,],z=Ztilde,index=2:41,marginal=FALSE,loo=T))
  ## -295.6863
  
  log.lik(sample.mcmc$all_K[1,,],z=Ztilde,index=2:41,marginal=TRUE)
  sum(log.lik(sample.mcmc$all_K[1,,],z=Ztilde,index=2:41,marginal=TRUE,loo=T))
  ## -301.1188
  
  log.lik(sample.mcmc$all_K[1,,],z=Ztilde,index=1,marginal=FALSE)
  sum(log.lik(sample.mcmc$all_K[1,,],z=Ztilde,index=1,marginal=FALSE,loo=T))
  ## - 6.298
}


## helper function to derive log-likelihood for each individual row of data
log.lik.loo <- function(theta,Ztilde,mu=NULL){
  ## theta: precision matrix
  ## mu : mean vector
  ## Ztilde: observed data as a matrix (even if with one column)
  ## value: log-likhood for each individual row of data
  prod <- rep(NA,nrow(Ztilde))
  
  ## SSE
  n <- nrow(Ztilde)
  p <- ncol(Ztilde)
  ##mu
  if(is.null(mu)){
    mu <- rep(0,p)
  }
  a2 <- Ztilde-mu
  #a3 <- apply(a2,1,function(x) outer(x,x,FUN="*"))
  #if(is.null(dim(a3))){
  #  S <- matrix(sum(a3),p,p)
  #}else{
  #  S <- matrix(apply(a3,1,sum),p,p)
  #}
  #prod__ <- crossprod(theta,S)
  #llik <- 0.5*n*log(det(theta)) - 0.5*sum(diag(prod__))
  #cat("log likelihood=",llik,"\n")
  
  logdet <- log(det(theta))
  for(i in 1:nrow(Ztilde)){
    prod[i] <- 0.5*(logdet - sum(a2[i,]*crossprod(theta,a2[i,])))
  }
  prod
}

