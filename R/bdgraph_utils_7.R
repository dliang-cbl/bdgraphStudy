#require(BDgraph)
## v. 7 handles predictive inference for near continuous non-Gaussian data
## 06-27-2018

## alternatively, always using single core, and use script language
## to run the multiple chains, this maybe easier.
## collect results from the disk then.
## 06-15-18: add thin argument to save only every other iterations
##           more-over store the thinned results there
## 06-20-18: add predictive distribution for model diagnostics and validation
## 05-06-2020: add a flag whether to predict
## value: bdgraph object plus all_Y: the predictive distribution
bdgraph.helper <- function(
  data, n=NULL, method="ggm",
  iter = 5000, burnin = iter/2, thin=1,
  g.start=NULL,g.space=NULL, g.prior=0.5,
  prior.df = 3, multi.update=NULL, print=1000, pp=TRUE,
  cores = 1
){

  ## rank data and dimension
  n <- dim(data)[1]
  p <- dim(data)[2]
  R <- 0 * as.matrix(data)
  for( j in 1:p ) R[ , j ] = match( data[ , j ], sort( unique( data[ , j ] ) ) ) 
  R[ is.na(R) ] = 0     # dealing with missing values	
  
  ## initial latent variable
  if((!any(is.na(data))) && (min( apply( R, 2, max ) ) > ( n - 5 * n / 100 ))){
    ## near continuous data
    Z <- qnorm( apply( data, 2, rank ) / ( n + 1 ) )
    is_near_continuous <- TRUE
  }else{
    # for non-Gaussian data
    is_near_continuous <- FALSE
    Z                  <- qnorm( apply( data, 2, rank, ties.method = "random" ) / ( n + 1 ) )
    Zfill              <- matrix( rnorm( n * p ), n, p )       # for missing values
    Z[ is.na( data ) ] <- Zfill[ is.na( data ) ]               # for missing values
    Z                  <- t( ( t( Z ) - apply( Z, 2, mean ) ) / apply( Z, 2, sd ) )
  }
  
  
  ## load shared library to sample from full conditional of latent variable
  #catch <- try(dyn.load("./BDgraph_src/check_os"),silent = TRUE)
  #if(class(catch)=="try-error"){
  #  dyn.load("./BDgraph_src/check_os.so")
  #}
  
  ## initial adjacency matrix
  ## sampler initialized by randomly sampling from the
  ## g-Wishart distribution with prior adjacency
  if(is.null(g.start)){
    g.start <- bdgraph.init(p=ncol(data),g.space = g.space,print=TRUE)
  }
  
  ## define batches
  sel <- seq(burnin+1,iter,thin) ## iterations to save
  iters__ <- diff(sel)  ## number of iterations for each batch
  
  ## make storage of all graphs and weights
  all_graphs_raw <- rep("",length(sel)) ## all graphs
  all_weights_raw <- rep(NA,length(sel)) ## corresponding weights 
                          ##(see version 5 discussion for details)
  
  ## run the burn-in
  cat("model is updating ",burnin, " burnin iterations.\n")
  
  currbatch__ <- bdgraph.study(
      data=data,n=n,method=method,algorithm="rjmcmc",iter=burnin+1,
      burnin=burnin,g.start=g.start,g.space=g.space,g.prior=g.prior,
      prior.df=prior.df,multi.update=multi.update,save.all=TRUE,
      print=print,cores=cores)
  
  ## determine whether data were treated as continuous
  is_gcgm__ <- ifelse(currbatch__$method=="gcgm",TRUE,FALSE)
  
  ## storage for K matrices
  ## the posterior sample of concentration matrix
  all_k <- array(NA,c(length(sel),dim(currbatch__$last_K))) 
  
  ## storage for latent variables
  ## the posterior samples of latent variable matrix (large!)
  all_z <- array(NA,c(length(sel),dim(currbatch__$last_Z)))
  

  ## store the first post-burnin iteartion
  all_graphs_raw[1] <- currbatch__$sample_graphs[1]
  all_weights_raw[1] <- currbatch__$graph_weights[1]
  all_k[1,,] <- currbatch__$last_K
  
  ## generate or retrieve latent variables
  if((!is_gcgm__) || any(is.na(currbatch__$last_Z))) {
    ## either ggm approximation used or NA in latent samples
    set.seed(round(100000*runif(1))) ## set random seed
    ## update Z from current concentration matrix estimate
    curr_Z__ <- .C("copula",Z=as.double(Z),K=as.double(currbatch__$last_K),
       R=as.integer(R),n=as.integer(n),p=as.integer(p),PACKAGE = "bdgraphStudy")$Z
    currbatch__$last_Z <- matrix(curr_Z__,n,p)
  }
  all_z[1,,] <- currbatch__$last_Z
  
  if(is_near_continuous){ 
    ## use data directly to be consistent with bdgraph
    currbatch__$last_Z <- Z
  }
  
  ## reporting module to make consistent with print
  cum_iters__ <- burnin+cumsum(iters__)  ## all iterations
  cum_print__ <- floor(cum_iters__/print) ## number of print batches
  to_print__ <- c(T,(diff(cum_print__)>0)) ## whether iterations go over
  print_iters__ <- cum_print__ * print ## iterations to print

  ## run the post burn in batches
  for(i in 1:length(iters__)){
    if(to_print__[i]){
      cat("iterations ",print_iters__[i],"\n")
    }
    ## restart the batch
    currbatch__ <- bdgraph.study(
      data=data,n=n,method=method,algorithm="rjmcmc",iter=iters__[i],
      burnin=0, g.start = currbatch__, g.space=g.space,g.prior=g.prior,
      prior.df=prior.df,multi.update=multi.update,save.all=TRUE,
      print=1e6,cores=cores)
    
    ## store each additional batch
    last_graph_ <- with(currbatch__, sample_graphs[all_graphs[length(all_graphs)]])
    all_graphs_raw[i+1] <- last_graph_
    all_weights_raw[i+1] <- with(currbatch__,all_weights[length(all_weights)])
    all_k[i+1,,] <- currbatch__$last_K
    if((!is_gcgm__)||any(is.na(currbatch__$last_Z))){
      seed <- 1000000*runif(1)
      set.seed(seed)
      curr_Z__ <- .C("copula",Z=as.double(all_z[i,,]),
                     K=as.double(currbatch__$last_K),
                     R=as.integer(R),n=as.integer(n),p=as.integer(p),
                     PACKAGE = "bdgraphStudy")$Z
      if(FALSE){
        rm(list=ls())
        load(file="Scratch/ggm_z6.rda")
        
        dyn.load("BDgraph_src/check_os.so")
        sink("BDgraph_src/ggm_debug_r6.txt")
        set.seed(seed)
        curr_Z__ <- .C("copula",Z=as.double(z),K=as.double(k),
                       R=as.integer(r),n=as.integer(n),p=as.integer(p),
                       PACKAGE = "bdgraphStudy")$Z
        dyn.unload("BDgraph_src/check_os.so")
        sink()
        range(curr_Z__)
      }
      if(any(is.na(curr_Z__))){
        save.image(file="tmp.rda")
        stop("missing values in latent variables.\n")
      }
      currbatch__$last_Z <- matrix(curr_Z__,n,p)
    }
    all_z[i+1,,] <- currbatch__$last_Z
    if(is_near_continuous){
      currbatch__$last_Z <- Z
    }
  }
  
  ## assemble the results
  sample_graphs <- unique(all_graphs_raw)
  all_graphs <- match(all_graphs_raw,sample_graphs)
  all_weights <- all_weights_raw
  graph_weights <- tapply(all_weights,all_graphs,sum)
  last_graph <- currbatch__$last_graph
  last_K <- currbatch__$last_K
  K_hat <- apply(all_k,c(2,3),mean) ## based on rj, see notes for details

  all_y <- NULL
  if(pp){
    cat("Updating done, predicting.\n")
    all_y <- all_z
    ## generate predictive distribution through inverse copula transfer
    for(m in 1:dim(all_z)[1]){
      sigma <- sqrt(diag(solve(all_k[m,,])))
      zs <- t(t(all_z[m,,])/sigma) ## standardize latent variables
      ps <- pnorm(zs)
      for(j in 1:ncol(zs)){
        all_y[m,,j] <- quantile(data[,j],probs=ps[,j],type=1,na.rm=TRUE)
      }
    }    
  }

  ## un-load shared library
  #if(class(catch)=="try-error"){
  #  dyn.unload("./BDgraph_src/check_os.so")
  #}
  #else{
  #  dyn.unload("./BDgraph_src/check_os")
  #}
  
  cat("Predicting done. \n")
  ##
  output <- list(sample_graphs=sample_graphs,all_graphs=all_graphs,
              all_weights=all_weights,graph_weights=graph_weights,
              K_hat=K_hat,all_K=all_k,all_Y=all_y,all_Z=all_z,Z=Z,
              method=ifelse(is_gcgm__,"gcgm","ggm"),
              last_graph=last_graph,last_K=last_K)
  
  class(output) <- "bdgraph"
  return(output)
}

test7 <- function()
{
  rm(list=ls())
  ## Session > Work Dir > Source File Loc
  library(BDgraph)
  ## Load un-exported R scripts from the package
  wd <- setwd("./BDgraph_R/")
  sink <- sapply(list.files(),source)
  setwd(wd)
  ## Load a local version of the main function (edited)
  source("bdgraph_study.R")
  
  ## simulate data
  set.seed(123)
  data.sim <- bdgraph.sim( n = 20, p = 6, size = 7, vis = TRUE )
  dat <- data.sim$data
  
  ## Load helper script to be debugged
  source("bdgraph_utils_7.R")
  sample <- bdgraph.helper(
    data = dat, method="gcgm",iter=150000,burnin=25000,
    thin=500,cores=1,print=1000)
  sample$method
  summary(sample)
  compare(data.sim,sample,colnames=c("True graph","BDgraph"))
  
  plot(sample$all_K[,1,2],type="l")
  acf(sample$all_K[,1,2])
  plot(sample$all_Y[,1,2],type="l")
  acf(sample$all_Y[,1,2])
  any(is.na(sample$all_Z))
  ## FALSE
  plot(dat[,1],sample$all_Y[1,,1],xlab="Observed",
       ylab="Posterior draw (iter=1)",
       main=colnames(dat)[1])
  ## all over the place!
  any(is.na(sample$all_Y))
  ## FALSE
  ## posterior median
  fity <- apply(sample$all_Y,c(2,3),median)
  op <- par(mfrow=c(3,2),mar=c(3.5,3.5,3.5,1.1))
  for(i in 1:6){
    plot(dat[,i],fity[,i],xlab="Observed",
         ylab="Posterior median",main=colnames(dat)[i])
    cor.test(fity[,i],dat[,i],method="kendall")
  }
  par(op)
  
}

test7b <- function(){
  rm(list=ls())
  ## Session > Work Dir > Source File Loc
  load(file="Scratch/tmp.rda")
  ## load the package
  library(BDgraph)
  ## Load un-exported R scripts from the package
  wd <- setwd("./BDgraph_R/")
  sink <- sapply(list.files(),source)
  setwd(wd)
  ## Load a local version of the main function (edited)
  source("bdgraph_study.R")
  ## Load helper script to be debugged
  source("bdgraph_utils_7.R")
  # set.seed(123456)
  # sample <- bdgraph.helper(
  #   data = dat, method="gcgm", g.space=g.space,
  #   iter=200000,burnin=0,thin=5,cores=1,print=1000)
  ## worked?
  set.seed(123456)
  sample <- bdgraph.helper(
    data = dat, method="gcgm", g.space=g.space,
    iter=2000000,burnin=0,thin=5000,cores=1,print=1000)
  
  sample$method
  plot(sample$all_K[,1,2],type="l")
  acf(sample$all_K[,1,2])
  plot(sample$all_Y[,1,2],type="l")
  acf(sample$all_Y[,1,2])
  any(is.na(sample$all_Z))
  plot(dat[,1],sample$all_Y[1,,1],xlab="Observed",
       ylab="Posterior draw (iter=1)",
       main=colnames(dat)[1])
  fity <- apply(sample$all_Y,c(2,3),median)
  op <- par(mfrow=c(3,3),mar=c(3.5,3.5,3.5,1.1))
  for(i in 1:9){
    plot(dat[,i],fity[,i],xlab="Observed",
         ylab="Posterior median",main=colnames(dat)[i])
    cor.test(fity[,i],dat[,i],method="kendall")
  }
  par(op)
  summary(sample)
  
  fit <- apply(sample$all_Z,c(2,3),quantile,probs=c(0.025,0.975))
  matplot(cbind(sample$Z[,1],fit[1,,1],fit[2,,1]),xlab="distance",
       ylab="copula transfer",main=colnames(dat)[1],type="l",
       lwd=2,lty=1,col=c(2,1,1))
}
## utility code to define a graphical space to explore
## loc: a 2 column matrix defining the unique links to explore
## p  : a dimension of the matrix
## value: a pxp adjacency matrix required by bdgraph code
bdgraph.space <- function(loc,p){
  stopifnot(max(c(loc))<=p)
  ## wish list remove duplicates
  raw <- matrix(0,p,p) ## raw matrix
  raw[loc] <- 2
  (raw + t(raw))/2
}

bdgraph.show <- function(sample.graph,dim){
  adj <- matrix(0,dim,dim)
  link <- as.numeric(strsplit(sample.graph,split="")[[1]])
  adj[upper.tri(adj)] <- link
  adj[lower.tri(adj)] <- t(adj)[lower.tri(adj)]
  adj
}
test <- function()
{
  diag(k) <- 0
  for(i in 1:41){
    cat("i=",i,"\n")
    print(range(k[i,][m[i,]==0]))
    print(range(abs(k[i,])[m[i,]==1]))
  }
}
## randomly sample from the graph space of size 2^(p*(p-1)/2)
## p: dimension of the graph
## g.space: sub-space of the graph to search, default is the whole graph
## value: pxp adjacency matrix denoting random links between nodes
bdgraph.init <- function(p,g.space=NULL,print=FALSE){
  ## all possible links
  num.pairs <- choose(p,2)
  ## generate random links 
  set.seed(as.numeric(Sys.time()))
  u <- runif(num.pairs)
  link <- as.numeric(u > 0.5)
  ## create storage matrix
  adj <- matrix(0,p,p)
  ## assign links
  adj[upper.tri(adj)] <- link
  adj[lower.tri(adj)] <- t(adj)[lower.tri(adj)]
  
  ## restrict to subspace
  if(!is.null(g.space)){
    adj[g.space<1] <- 0
    if(print){
      graph <- paste(g.space[upper.tri(g.space)],collapse = "")
      cat("restrict space=",graph,"\n")
    }
  }
  
  if(print){
    graph <- paste(adj[upper.tri(adj)],collapse = "")
    cat("initial graph =",graph,"\n")
  }
  ## return
  adj
}


require(coda)
## helper function to load results from multiple bdgraph objects
## x: a list of "rda" file with stored sample
## burnin: additional burnin to dischard, default zero
## y: character string the name of the stored object
## return:
##  value : the concatenated bdgraph object
##  mcmc  : a mcmc object for the size of the graph
##
## see bdgraph_utils.docx
bdgraph.load <- function(x,burnin=0,thin=1,y="sample")
{
  results.files <- x
  
  ## initialize dimensions
  load(file=results.files[1])
  sample__ <- get(y)
  p <- ncol(sample__$last_graph) ## graph dimension
  qp <- choose(p,2)  ## all links
  M <- length(results.files)   ## number of parallel chains
  niter <- length(sample__$all_graphs) ## number of iterations
  sel <- seq(burnin+1,niter,thin)
  
  ## results storage
  sample.graphs.list <- vector("list",M)
  graph.weights.list <- vector("list",M)
  K.hat.list <- array(NA,c(M,p,p))
  all.graphs.list <- vector("list",M)
  all.weights.list <- vector("list",M)
  last.graph.list <- array(NA,c(M,p,p))
  last.K.list <- array(NA,c(M,p,p))
  dname <- c(list(NULL),dimnames(sample__$K_hat))
  dimnames(K.hat.list) <- dname
  dimnames(last.graph.list) <- dname
  dimnames(last.K.list) <- dname
  
  all.K.list <- vector("list",M)
  all.Y.list <- vector("list",M)
  all.Z.list <- vector("list",M)
  
  ## convergency check storage
  size.list <- vector("list",M)
  expected.size.list <- vector("list",M)
  
  ## start collecting results
  for(i in 1:M){
    cat(results.files[i],"\n")
    ## Load stored object
    load(file=results.files[i])  
    sample__ <- get(y)
    
    ## store Monte Carlo Samples
    all.K.list[[i]] <- sample__$all_K[sel,,]
    all.Y.list[[i]] <- sample__$all_Y[sel,,]
    all.Z.list[[i]] <- sample__$all_Z[sel,,]
    
    ##browser()
    ## apply the burn_in -  weights aggregated
    all_graphs_raw__ <- with(sample__,sample_graphs[all_graphs[sel]])
    all_weights__ <- sample__$all_weights[sel]
    sample_graphs__ <- unique(all_graphs_raw__)
    all_graphs__ <- match(all_graphs_raw__,sample_graphs__)
    #graph_weights__ <- as.numeric(table(all_graphs__))
    graph_weights__ <- tapply(all_weights__,all_graphs__,sum)
    
    
    ## calculate graph size as a check for convergence
    #bdgraph.obj <- sample__
    #all_graphs = bdgraph.obj$all_graphs[sel]
    #sample_graphs <- bdgraph.obj$sample_graphs
    sizesample_graphs <- sapply(sample_graphs__,function(x) length(which(unlist(strsplit(as.character(x),""))==1)))
    sizeall_graphs = sizesample_graphs[all_graphs__]
    ## store the size
    size.list[[i]] <- sizeall_graphs
    
    #browser()
    ## what is the posterior expected number of edges at each iteration?
    ## calculate expected size as the sum of posterior probability
    ## of all links, where posterior probability of each link
    ## is the cumulative proportion of included links
    dump <- matrix(0,qp,length(all_graphs_raw__))
    vec_result <- rep(0,qp)
    for(g in 1:length(all_graphs_raw__)){
      which.edge <- which(unlist(strsplit(as.character(all_graphs_raw__[g]),"")) == 1)
      vec_result[which.edge] <- vec_result[which.edge]+all_weights__[g]
      dump[,g] <- vec_result/sum(all_weights__[c(1:g)])
    }
    expected.size.list[[i]] <- apply(dump,2,sum)
    
    ## assuming convergence collect the results
    sample.graphs.list[[i]] <- sample_graphs__
    graph.weights.list[[i]] <- graph_weights__
    K.hat.list[i,,] <- sample__$K_hat
    all.graphs.list[[i]] <- all_graphs__
    all.weights.list[[i]] <- all_weights__
    last.graph.list[i,,] <- sample__$last_graph
    last.K.list[i,,] <- sample__$last_K
  }
  
  ## make convergence output
  tmplist <- vector("list",M)
  for(i in 1:M){
    tmplist[[i]] <- as.mcmc(size.list[[i]])
  }
  size.mcmc <- as.mcmc.list(tmplist)
  
  tmplist <- vector("list",M)
  for(i in 1:M){
    tmplist[[i]] <- as.mcmc(expected.size.list[[i]])
  }
  expected.size.mcmc <- as.mcmc.list(tmplist)
  
  if(F){
    sample.graphs.list <- list(
      c("a","b","c"),
      c("a","b","d")
    )
    all.graphs.list <- list(
      c(1,2,1,1,3),
      c(1,1,1,2,3)
    )
    M <- 2
  }
  ## make sample graph object
  sample_graphs_all <- do.call(c,sample.graphs.list)
  sample_graphs__ <- unique(sample_graphs_all) ## unique sample_graphs
  
  ## make all graph object
  tmp.allgraph <- vector("list",M)
  for(i in 1:M){
    #browser()
    tmp_sample_graph_ <- match(sample.graphs.list[[i]],sample_graphs__)
    tmp.allgraph[[i]] <- tmp_sample_graph_[all.graphs.list[[i]]]
  }
  all_graphs__ <- do.call(c,tmp.allgraph)
  
  if(FALSE){
    ## test code
    g0.list <- vector("list",M)
    for(i in 1:M){
      g0.list[[i]] <- sample.graphs.list[[i]][all.graphs.list[[i]]]
    }
    g0 <- do.call(c,g0.list)
    g1 <- sample_graphs__[all_graphs__]
    identical(g0,g1)
    ## TRUE debug issues
  }
  #browser()
  ## make graph weights objects - assuming weighted
  #graph_weights__ <- as.numeric(table(all_graphs__))
  graph_weights__ <- tapply(do.call(c,all.weights.list),all_graphs__,sum)
  

  ## make all weights objects for reversible jump MCMC
  #all_weights__ <- rep(1,length(all_graphs__))
  all_weights__ <- do.call(c,all.weights.list)
  
  ## make estimates
  K_hat__ <- apply(K.hat.list,c(2,3),mean)
  last_graph__ <- last.graph.list[M,,]
  last_K__ <- last.K.list[M,,]
  
  dimnames(last_graph__) <- dimnames(sample__$last_graph)
  dimnames(last_K__) <- dimnames(sample__$last_K)
  
  ## make posterior summaries
  bdgraph.obj <- list(
    sample_graphs=sample_graphs__,
    graph_weights=graph_weights__,
    K_hat=K_hat__,
    all_graphs=all_graphs__,
    all_weights=all_weights__,
    last_graph=last_graph__,
    last_K = last_K__,
    Z = sample__$Z,
    all_K = abind(all.K.list,along=1),
    all_Z = abind(all.Z.list,along=1),
    all_Y = abind(all.Y.list,along=1),
    size_mcmc=size.mcmc,
    emcmc=expected.size.mcmc
  )
  class(bdgraph.obj) = "bdgraph"
  
  #list(value=bdgraph.obj,mcmc=size.mcmc,emcmc=expected.size.mcmc)
  bdgraph.obj
}

test3 <- function()
{
  rm(list=ls())
  ## on cluster Fall2017/BIGDATAJellyfish/cggm7
  results.files <- list.files(path="./cggm7c",pattern = "res",full=TRUE)
  source("bdgraph_utils_4.R")
  #load(file="./cggm7c/res_2_1.rda")
  fn <- c("./cggm7c/res_2_1.rda","./cggm7c/res_2_2.rda")
  test.results <- bdgraph.load(x=fn,y="sample",thin=5000)
  summary(test.results$value)
}

## x: fitted bdgraph object
## name: variable names (default to built in)
## th: threshold to determine posterior link
##     if null, try to use the most likely graph
cggm.bdgraph.Cgraph <- function(x,name=NULL,th=0.5)
{
  if(is.null(th)){
    results__ <- pgraph(x,number.g=1)
  }else{
    results__ <- select(x,cut=th,vis=FALSE)
  }
  
  results__ <- as.matrix(results__)
  ## find the edge selected
  sel <- which(results__[upper.tri(results__)]>0)
  node1 <- row(results__)[upper.tri(results__)][sel]
  node2 <- col(results__)[upper.tri(results__)][sel]
  ## find the fitted concentration matrix
  fittedConc <- x$K_hat[upper.tri(results__)][sel]
  ## check results
  stopifnot(all(!is.na(sel)))
  stopifnot(all(!is.na(fittedConc)))
  ## stroage of results
  mat <- matrix(0,nrow(results__),ncol(results__))
  if(is.null(name)){
    dimnames(mat) <- dimnames(x$K_hat)
  }
  else{
    dimnames(mat) <- list(name,name)
  }
  if(length(sel)>0){
    val <- ifelse(fittedConc>0,100,10)
    mat[cbind(node1,node2)] <- val
    mat[cbind(node2,node1)] <- val
  }
  mat
}


require(tidyr)
require(mgcv)
## helper function to translate a gam formula into a regular lm formula
## x: a formula possibly with s terms on the RHS
## value: a formula with s terms replaced with linear term
bdgraph.formula <- function(x){
  l.terms <- as.character(x)
  l.loc <- switch(length(l.terms)-1,2,3)
  l.resp <- switch(length(l.terms)-1,NULL,l.terms[2])
  l.formula <- reformulate(l.terms[l.loc])
  l.pred <- all.vars(l.formula)
  if(length(l.pred)==0){
    ## intercept only
    l.pred <- "1"
  }
  as.formula(paste(l.resp,"~",l.pred,sep=""))
}
test <- function()
{
  source("bdgraph_utils_2.R")
  bdgraph.formula(cbind(y,z)~s(x))
  bdgraph.formula(cbind(y,z)~x+w)
  bdgraph.formula(~x+w)
  bdgraph.formula(y~1)
  bdgraph.formula(cbind(y,z)~1)
}
## spread does require an id variable, 
bdgraph.detrend <- function(formula,data,bytime=FALSE){
  mf <- model.frame(bdgraph.formula(formula),data=data) 
  ncruise <- ncol(mf[[1]]) ## default response number of cruises
  ## gather cruise data
  tmp <- gather(cbind(data[all.vars(formula)],rid=1:nrow(data)),"cruise__","y__",1:ncruise)
  ## develop trend formula
  if(bytime){
    trend <- reformulate(c("cruise__",as.character(formula)[3]),response = "y__")
  }else{
    trend <- reformulate(as.character(formula)[3],response = "y__")
  }
  ## detrending
  trend.gam <- gam(trend,data=tmp)
  ## residuals
  tmp[,"y__"] <- resid(trend.gam)
  ## spread the response out again
  out <- spread(tmp,cruise__,y__)
  #browser()
  ## re-order the rows
  oo <- order(out[,"rid"])
  ## re-order the columns
  r <- out[oo,all.vars(formula)]
  
  list(value=r,raw=tmp,fitted=trend.gam)
}

test4 <- function()
{
  rm(list=ls())
  gc(T)
  ## load spatial aligned data
  load(file="../Data/paxsur_4f.rda")
  source("bdgraph_utils_2.R")
  creek <- subset(paxsur,transect=="SLCV")
  fout <- bdgraph.detrend(formula=cbind(N.4,N.5,N.6,N.7)~s(cmin_d),
                          data=creek)
  plot(creek$cmin_d,creek$N.5)
  plot(creek$cmin_d,fout$value$N.5)
}