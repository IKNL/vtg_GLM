node_beta=function(formula,Data,weights=NULL,family=gaussian,beta=NULL,iter){
  y=eval(formula[[2]], envir = Data) #extract y and X varibales name from formula
  X=model.matrix(formula,data = Data) #create a model matrix 
  offset=model.offset(model.frame(formula,data = Data)) #extract the offset from formula (if exists)
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
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

master_beta=function(...,beta=NULL,family=gaussian){
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
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
  if(family$family %in% c('poisson','binomial')){
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

node_deviance=function(formula,Data,weights=NULL,family=gaussian,beta,iter){
  y=eval(formula[[2]], envir = Data)
  X=model.matrix(formula,data = Data)
  offset=model.offset(model.frame(formula,data = Data)) 
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
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

master_convergence=function(...,tol= 1e-08,beta,iter,family=gaussian,formula,maxit=25){
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  x <- list(...)
  dev_old=Reduce(`+`,lapply(1:length(x),function(j) x[[j]]$dev_old))
  dev=Reduce(`+`,lapply(1:length(x),function(j) x[[j]]$dev))
  convergence=(abs(dev - dev_old) / (0.1 + abs(dev)) < tol)
  if(convergence==FALSE & iter<maxit){
    return(list(converged=convergence))
  }else{
   zvalue <- beta$coef[,ncol(beta$coef)]/beta$se
    if(beta$est.disp){
      pvalue <- 2 * pt(-abs(zvalue), beta$nobs-beta$nvars)
    }else{
      pvalue <- 2 * pnorm(-abs(zvalue))
    }
   L=list(converged=convergence,
          coefficients=beta$coef[,ncol(beta$coef)],
          Std.Error=beta$se,
          pvalue=pvalue,
          zvalue=zvalue,
          dispersion=beta$disp,
          est.disp=beta$est.disp,
          formula=formula,
          family=family,
          iter=iter,
          deviance=dev,
          null.deviance=dev,
          nobs=beta$nobs,
          nvars=beta$nvars)
  return(L)
  }
 }

as.GLM <- function(obj, data=NULL) {
  
  dots <- as.list(obj$coefficients)
  
  out<-list()
  
  tt <- terms(obj$formula, data=data)
  
  if(!is.null(data)) {
    
    mf <- model.frame(tt, data)
    
    vn <- sapply(attr(tt, "variables")[-1], deparse)
    
    
    tt=M_2$terms
    vn <- sapply(attr(tt, "variables")[-1], deparse)
    
    
    
    if((yvar <- attr(tt, "response"))>0)
      
      vn <- vn[-yvar]
    
    xlvl <- lapply(data[vn], function(x) if (is.factor(x))
      
      levels(x)
      
      else if (is.character(x))
        
        levels(as.factor(x))
      
      else
        
        NULL)
    
    attr(out, "xlevels") <- xlvl[!vapply(xlvl,is.null,NA)]
    
    attr(tt, "dataClasses") <- sapply(data[vn], stats:::.MFclass)
    
  }
  
  out$terms <- tt
  
  coef <- numeric(0)
  
  stopifnot(length(dots)>1 & !is.null(names(dots)))
  
  for(i in seq_along(dots)) {
    
    if((n<-names(dots)[i]) != "") {
      
      v <- dots[[i]]
      
      if(!is.null(names(v))) {
        
        coef[paste0(n, names(v))] <- v
        
      } else {
        
        stopifnot(length(v)==1)
        
        coef[n] <- v
        
      }
      
    } else {
      
      coef["(Intercept)"] <- dots[[i]]
      
    }   
    
  }
  
  out$coefficients <- coef
  
  out$rank <- length(coef)
  
  family=obj$family
  
  if (!missing(family)) {
    
    out$family <- if (class(family) == "family") {
      
      family
      
    } else if (class(family) == "function") {
      
      family()
      
    } else if (class(family) == "character") {
      
      get(family)()
      
    } else {
      
      stop(paste("invalid family class:", class(family)))
      
    }
    
    out$qr <- list(pivot=seq_len(out$rank))
    
    out$deviance <- obj$deviance
    
    out$null.deviance <- obj$null.deviance
    
    out$aic <- 1
    
    out$iter=obj$iter
    
    out$df.null=obj$nobs-1
    
    out$df.residual=obj$nobs-out$rank
    
    out$call=call("glm_FL",f,family=obj$family$family)
    
    class(out) <- c("glm","lm")
    
    
  } else {
    
    class(out) <- "lm"
    
    out$fitted.values <- predict(out, newdata=dd)
    
    out$residuals <- out$mf[attr(tt, "response")] - out$fitted.values
    
    out$df.residual <- nrow(data) - out$rank
    
    out$model <- data
    
    #QR doesn't work
    
  }
  
  out
  
}

Summary_FL_GLM <- function(obj, ...){
  sign=cut(obj$pvalue,c(0,0.001,0.01,0.05,0.1,1),include.lowest = T,labels=c('***','**','*','.',''))
    cat("Model:\n")
  print(obj$formula)
  
  DF=data.frame(round(obj$coefficients,5),round(obj$Std.Error,5),round(obj$zvalue,3),round(obj$pvalue,5),sign)
  if(obj$est.disp){
    names(DF)=c("Estimate","Std. Error","t value","Pr(>|t|)","")
  }else{
    names(DF)=c("Estimate","Std. Error","z value","Pr(>|z|)","")
  }
  cat("\nCoefficients:\n")
  print(DF)
  
  cat("---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n\n")
  cat(paste0('(Dispersion parameter for ',obj$regression," family taken to be ",round(obj$dispersion,5),')\n \n'))
  cat(paste0("Residual deviance: ",round(obj$deviance,2),"  on ",obj$nobs-obj$nvars," degrees of freedom \n \n"))
  cat(paste0('Number of Fisher Scoring iterations: ',obj$iter,'\n'))

}


# 
# glm_obj=list(
#   call=paste0('glm_FL(formula = ',Reduce(paste, deparse(formula)),', family = ',reg_type,')'),
#   coefficients=beta$coef[,ncol(beta$coef)],
#   rank=beta$nvars+1,
#   family=reg_type,
#   deviance=dev,
#   null.deviance=dev,
#   iter=iter,
#   df.null=beta$nobs-1,
#   df.residual=beta$nobs-beta$nvars,
#   converged=convergence,
#   formula=formula,
#   aic=0,
#   control=list(epsilon=tol,maxit=maxit),
#   terms=terms.formula(f)
# ) 
# L=structure(glm_obj, class = "glm")


# 
# 
