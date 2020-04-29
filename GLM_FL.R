master_init=function(formula,family=gaussian,tol= 1e-08,maxit=25){
  saveRDS(list(formula=formula,family=family,iter=1,tol= tol,maxit=maxit),"master.Rds")
}

node_beta=function(Data,weights=NULL,master=NULL,userID){
  #if(is.null(user)){ print("Please specify the user number (1,2,3,....)"); break}
  formula=master$formula 
  family=master$family 
  y=eval(formula[[2]], envir = Data) #extract y and X varibales name from formula
  X=model.matrix(formula,data = Data) #create a model matrix 
  offset=model.offset(model.frame(formula,data = Data)) #extract the offset from formula (if exists)
  #functions of the family required (gaussian, poisson, logistic,...)
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
  
  if(master$iter==1){
    etastart = NULL
    eval(family$initialize) # initializes n and fitted values mustart
    eta = family$linkfun(mustart) # we then initialize eta with this
  } else {
    eta = (X %*% master$coef[,ncol(master$coef)]) + offset #update eta
  }

  mu = family$linkinv(eta) 
  varg = family$variance(mu)
  gprime = family$mu.eta(eta)
  z = (eta - offset) + (y - mu) / gprime #calculate z
  W = weights * as.vector(gprime^2 / varg) #update the weights
  dispersion <- sum(W *((y - mu)/family$mu.eta(eta))^2) #calculate the dispersion matrix
  output=list(v1=crossprod(X,W*X),v2=crossprod(X,W*z),dispersion=dispersion,
              nobs=nobs,nvars=nvars,wt1=sum(weights * y),wt2=sum(weights))
  saveRDS(output,file = paste0(userID,'.Rds'))
  #return(output)
}

master_beta=function(...,nodes=NULL,master=NULL){
  #receive as many object as many are the nodes involved in the analysis (...)
  #the function update the betas
  formula=master$formula 
  family=master$family 
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if(is.null(nodes))  {
    g <- list(...)  #place the dots into a list
  }else{
    g <- nodes   
  }
  allwt=Reduce(`+`,lapply(1:length(g),function(j) g[[j]]$wt2)) #total sum o fweights
  wtdmu=Reduce(`+`,lapply(1:length(g),function(j) g[[j]]$wt1/allwt)) #global weighted mu
  a=Reduce(`+`,lapply(1:length(g),function(j) g[[j]]$v1)) #sum up components of the matrix to be inverted calculated in each node
  b=Reduce(`+`,lapply(1:length(g),function(j) g[[j]]$v2)) #sum up components of the matrix to be inverted calculated in each node
  phi=Reduce(`+`,lapply(1:length(g),function(j) g[[j]]$dispersion)) # sum up components dispersion matrix
  nobs=Reduce(`+`,lapply(1:length(g),function(j) g[[j]]$nobs)) #total number of observation
  nvars=nrow(g[[1]]$v1) #number of variables

  if(is.null(master)){
    beta = rep(1, nvars)
  }else{
    beta=master$coef
  }
  if(family$family %in% c('poisson','binomial')){
      disp=1
      est.disp=FALSE
    }else{
      disp=phi/(nobs-nvars)
      est.disp=T
    }
 
  fb=solve(a,b, tol=2*.Machine$double.eps) #calculate the new betas
  se=sqrt(diag(solve(a)*disp)) #calculate the Standard error of coefficients
 
    master$coef=cbind(master$coef,fb)
    master$se=se
    master$disp=disp
    master$est.disp=est.disp
    master$nobs=nobs
    master$nvars=nvars
    master$wtdmu=wtdmu
 
  saveRDS(master,file = paste0("master.Rds"))
  #return(output)
}

node_deviance=function(Data,weights=NULL,master,userID){
  #the function update the betas
  formula=master$formula 
  family=master$family 
  #the function calculate the residual deviance with updated betas for the single node
  y=eval(formula[[2]], envir = Data) #extract y variable names
  X=model.matrix(formula,data = Data) #extract X variables
  offset=model.offset(model.frame(formula,data = Data)) #extract the offset
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
  
  if(master$iter==1){ #only for first iteration
    etastart = NULL
    nobs = nrow(X)    # needed by the initialize expression below
    nvars = ncol(X)   # needed by the initialize expression below
    eval(family$initialize) # initializes n and fitted values mustart
    eta = family$linkfun(mustart) + offset # we then initialize eta
    mu_old = family$linkinv(eta)
    dev_old=0
  }else{
    mu_old = family$linkinv(X %*% master$coef[,ncol(master$coef)-1]) 
    dev_old=sum(family$dev.resids(y, mu_old,weights))
  }
  eta = X %*% master$coef[,ncol(master$coef)] +offset #calcaute updated eta
  mu = family$linkinv(eta-offset)  
  dev = sum(family$dev.resids(y, mu,weights)) #calculate new deviance
  dev.null=sum(family$dev.resids(y, master$wtdmu,weights))
  output=list(dev_old=dev_old,dev=dev,dev.null=dev.null)
  saveRDS(output,file = paste0(userID,".Rds"))
  #return(output)
}

master_deviance=function(...,nodes=NULL,master){
  #receive as many object as many are the nodes involved in the analysis (...)
  #the function evaluate if the algorithm converge
  #the function update the betas
  tol=master$tol
  maxit=master$maxit
  formula=master$formula
  family=master$family 
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if(is.null(nodes))  {
    x <- list(...)  #place the dots into a list
  }else{
    x <- nodes   
  }
  dev_old=Reduce(`+`,lapply(1:length(x),function(j) x[[j]]$dev_old)) #sum up deviance of previous iteration
  dev=Reduce(`+`,lapply(1:length(x),function(j) x[[j]]$dev)) #sum up new deviance
  dev.null=Reduce(`+`,lapply(1:length(x),function(j) x[[j]]$dev.null)) #sum up null deviance
  convergence=(abs(dev - dev_old) / (0.1 + abs(dev)) < tol) #evaluate if algorithm  converge
  if(convergence==FALSE & master$iter<maxit){
      master$converged=convergence
      master$iter=master$iter+1
 
    saveRDS(master,file = paste0("master.Rds"))
    #return(output)
  }else{
   zvalue <- master$coef[,ncol(master$coef)]/master$se
    if(master$est.disp){
      pvalue <- 2 * pt(-abs(zvalue), master$nobs-master$nvars)
    }else{
      pvalue <- 2 * pnorm(-abs(zvalue))
    }
   master=list(converged=TRUE,
          coefficients=master$coef[,ncol(master$coef)],
          Std.Error=master$se,
          pvalue=pvalue,
          zvalue=zvalue,
          dispersion=master$disp,
          est.disp=master$est.disp,
          formula=master$formula,
          family=family,
          iter=master$iter,
          deviance=dev,
          null.deviance=dev.null,
          nobs=master$nobs,
          nvars=master$nvars)
   saveRDS(master,file = paste0("master.Rds"))
   #return(output)
  }
  #return a list of output only if the algorithm converge
 }

as.GLM <- function(obj, data=NULL) {
  #fill a GLM object with output of the Federated Learning GLM
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

summary_FL_GLM <- function(obj, ...){
  #Summarize the output of the GLM in a user-friendly tableS
  sign=cut(obj$pvalue,c(0,0.001,0.01,0.05,0.1,1),include.lowest = T,labels=c('***','**','*','.',''))
    cat("Call:\n")
    print(call("glm_FL",obj$formula,family=obj$family$family))
 # print(obj$formula)
  
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
  cat(paste0('(Dispersion parameter for ',obj$family$family," family taken to be ",round(obj$dispersion,5),')\n \n'))
  cat(paste0("Null deviance: ",round(obj$null.deviance,2),"  on ",obj$nobs-1," degrees of freedom \n"))
  cat(paste0("Residual deviance: ",round(obj$deviance,2),"  on ",obj$nobs-obj$nvars," degrees of freedom \n \n"))
  cat(paste0('Number of Fisher Scoring iterations: ',obj$iter,'\n'))

}

FL_GLM=function(...,f,family,maxit=25,tol=1e-08){
  data=list(...)
  Master_1=NULL
  for(j in 1:maxit){
    Node_1=lapply(data,function(p) node_beta(formula = f,Data = p,beta = Master_1,iter = j,family =family))
    Master_1=master_beta(nodes=Node_1,beta = Master_1,family = family)
    Node_2=lapply(data,function(p) node_deviance(formula = f,Data = p,beta = Master_1,iter = j,family =family))
    Master_2=master_deviance(nodes=Node_2,beta=Master_1,iter=j,family = family,formula = f,maxit=maxit)
    if(Master_2$converged)break
  }
  return(Master_2)
}



                