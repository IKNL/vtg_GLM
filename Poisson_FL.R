node_beta=function(formula,Data,weights=NULL,reg_type=c('gaussian','poisson','binomial'),beta=NULL,iter){
  y=eval(formula[[2]], envir = Data)
  X=model.matrix(formula,data = Data)
  offset=model.offset(model.frame(formula,data = Data))
  #reg_type <- match.arg(reg_type)
  if(reg_type=='gaussian') family = gaussian()
  if(reg_type=='poisson')  family = poisson(log)
  if(reg_type=='binomial') family = binomial()
  if (is.null(weights)) weights <- rep.int(1, nrow(X))
  if (is.null(offset)) offset <- rep.int(0, nrow(X))
  nobs = nrow(X)    # needed by the initialize expression below
  nvars = ncol(X)   # needed by the initialize expression below
  if(iter==1){
    etastart = NULL
    eval(family$initialize) # initializes n and fitted values mustart
    eta = family$linkfun(mustart) # we then initialize eta with this
  } else {
    eta = (X %*% beta[,ncol(beta)]) + offset # potentially +offset if you would have an offset argument as well
  }
  mu = family$linkinv(eta) 
  varg = family$variance(mu)
  gprime = family$mu.eta(eta)
  z = (eta - offset) + (y - mu) / gprime # potentially -offset if you would have an offset argument as well
  W = weights * as.vector(gprime^2 / varg)
  dispersion <- sum(W *((y - mu)/family$mu.eta(eta))^2)
  
  return(list(v1=crossprod(X,W*X),v2=crossprod(X,W*z),dispersion=dispersion,nobs=nobs,nvars=nvars))
}

master_beta=function(...,beta=NULL,reg_type=c('gaussian','poisson','binomial')){
  g <- list(...)
  a=Reduce(`+`,lapply(1:length(g),function(j) g[[j]]$v1))
  b=Reduce(`+`,lapply(1:length(g),function(j) g[[j]]$v2))
  phi=Reduce(`+`,lapply(1:length(g),function(j) g[[j]]$dispersion))
  nobs=Reduce(`+`,lapply(1:length(g),function(j) g[[j]]$nobs))
  nvars=nrow(g[[1]]$v1)
  if(is.null(beta)){
    beta = rep(1, nvars)
  }else{
    beta=beta$coef
  }
  if(reg_type %in% c('poisson','binomial')){
      disp=1
      est.disp=FALSE
    }else{
      disp=phi/(nobs-nvars)
      est.disp=T
    }
  f=solve(a,b, tol=2*.Machine$double.eps)
  se=sqrt(diag(solve(a)*disp))
  return(list(coef=cbind(beta,f),se=se,disp=disp,est.disp=est.disp,nobs=nobs,nvars=nvars))
}

node_deviance=function(formula,Data,weights=NULL,reg_type=c('gaussian','poisson','binomial'),beta,iter){
  y=eval(formula[[2]], envir = Data)
  X=model.matrix(formula,data = Data)
  offset=model.offset(model.frame(formula,data = Data)) 
  if(reg_type=='gaussian') family = gaussian()
  if(reg_type=='poisson')  family = poisson(log)
  if(reg_type=='binomial') family = binomial()
  if (is.null(weights)) weights <- rep.int(1, nrow(X))
  if (is.null(offset)) offset <- rep.int(0, nrow(X))
  
  if(iter==1){
    etastart = NULL
    nobs = nrow(X)    # needed by the initialize expression below
    nvars = ncol(X)   # needed by the initialize expression below
    eval(family$initialize) # initializes n and fitted values mustart
    eta = family$linkfun(mustart) + offset# we then initialize eta with this
    mu_old = family$linkinv(eta)
    dev_old=0
  }else{
    mu_old = family$linkinv(X %*% beta[,ncol(beta)-1]) 
    dev_old=sum(family$dev.resids(y, mu_old,weights))
  }
  eta = X %*% beta[,ncol(beta)] +offset # potentially +offset if you would have an offset argument as well
  mu = family$linkinv(eta-offset) 
  dev = sum(family$dev.resids(y, mu,weights))
  return(list(dev_old=dev_old,dev=dev))
}

master_convergence=function(...,tol=1e-16,beta,iter,reg_type=c('gaussian','poisson','binomial'),formula,maxit){
  x <- list(...)
  dev_old=Reduce(`+`,lapply(1:length(x),function(j) x[[j]]$dev_old))
  dev=Reduce(`+`,lapply(1:length(x),function(j) x[[j]]$dev))
  convergence=(abs(dev - dev_old) / (0.1 + abs(dev)) < tol)
  if(convergence==FALSE & iter<maxit){
    return(list(convergence=convergence))
  }else{
   zvalue <- beta$coef[,ncol(beta$coef)]/beta$se
    if(beta$est.disp){
      pvalue <- 2 * pt(-abs(zvalue), beta$nobs-beta$nvars)
    }else{
      pvalue <- 2 * pnorm(-abs(zvalue))
    }
    L=list(convergence=convergence,
    coefficient=beta$coef[,ncol(beta$coef)],
    Std.Error=beta$se,
    pvalue=pvalue,
    zvalue=zvalue,
    dispersion=beta$disp,
    est.disp=beta$est.disp,
    model=formula,
    regression=reg_type,
    iter=iter,
    dev_res=dev,
    nobs=beta$nobs,
    nvars=beta$nvars)
  return(L)
  }
 }

Summary_FL_GLM <- function(x, ...){
  sign=cut(x$pvalue,c(0,0.001,0.01,0.05,0.1,1),include.lowest = T,labels=c('***','**','*','.',''))
    cat("Model:\n")
  print(x$model)
  
  DF=data.frame(round(x$coefficient,5),round(x$Std.Error,5),round(x$zvalue,3),round(x$pvalue,5),sign)
  if(x$est.disp){
    names(DF)=c("Estimate","Std. Error","t value","Pr(>|t|)","")
  }else{
    names(DF)=c("Estimate","Std. Error","z value","Pr(>|z|)","")
  }
  cat("\nCoefficients:\n")
  print(DF)
  
  cat("---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n\n")
  cat(paste0('(Dispersion parameter for ',x$regression," family taken to be ",round(x$dispersion,5),')\n \n'))
  cat(paste0("Residual deviance: ",round(x$dev_res,2),"  on ",x$nobs-x$nvars," degrees of freedom \n \n"))
  cat(paste0('Number of Fisher Scoring iterations: ',x$iter,'\n'))

}


