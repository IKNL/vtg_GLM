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
  return(list(coef=cbind(beta,f),se=se,disp=disp,est.disp=est.disp,rank=nobs-nvars))
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
  mu = family$linkinv(X %*% beta[,ncol(beta)]) 
  dev = sum(family$dev.resids(y, mu,weights))
  devnull = sum(family$dev.resids(y, sum(weights * y)/sum(weights),weights))
  return(list(dev_old=dev_old,dev=dev,devnull=devnull))
}

master_convergence=function(...,tol=1e-16,beta,iter,reg_type=c('gaussian','poisson','binomial'),formula){
  x <- list(...)
  dev_old=Reduce(`+`,lapply(1:length(x),function(j) x[[j]]$dev_old))
  dev_res=Reduce(`+`,lapply(1:length(x),function(j) x[[j]]$dev))
  dev_null=Reduce(`+`,lapply(1:length(x),function(j) x[[j]]$devnull))
  convergence=(abs(dev - dev_old) / (0.1 + abs(dev)) < tol)
  if(convergence==FALSE){
    return(list(convergence=convergence))
  }else{
   zvalue <- beta$coef[,ncol(beta$coef)]/beta$se
    if(beta$est.disp){
      pvalue <- 2 * pt(-abs(zvalue), beta$rank)
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
    dev_res=dev_res,
    dev_null=dev_null)
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
  cat(paste0('Number of Fisher Scoring iterations: ',x$iter,'\n'))
}


#sign=cut(pvalue,c(0,0.001,0.01,0.05,0.1,1),include.lowest = T,labels=c('***','**','*','.',''))
#cbind("Estimate"=beta$coef[,ncol(beta$coef)],"Std. Error"=beta$se,"z value"=zvalue,"Pr(>|z|)"=pvalue,""=sign)

###poisson example

p <- read.csv("https://stats.idre.ucla.edu/stat/data/poisson_sim.csv")
p <- within(p, {
  prog <- factor(prog, levels=1:3, labels=c("General", "Academic", 
                                            "Vocational"))
  id <- factor(id)
})

p1=p[1:(nrow(p)/2),-1]
p2=p[(nrow(p)/2+1):nrow(p),-1]

maxit=25 
tol=1e-16

M_1=NULL
f=num_awards ~ prog + math

#iterations
for(j in 1:maxit){
  N_1=node_beta(formula = f,Data = p1,beta = M_1$coef,iter = j,reg_type = 'poisson')
  N_2=node_beta(formula = f,Data = p2,beta = M_1$coef,iter = j,reg_type = 'poisson')
  M_1=master_beta(N_1,N_2,beta = M_1,reg_type = 'poisson')
  D_1=node_deviance(formula = f,Data = p1,beta = M_1$coef,iter = j,reg_type = 'poisson')
  D_2=node_deviance(formula = f,Data = p2,beta = M_1$coef,iter = j,reg_type = 'poisson')
  M_2=master_convergence(D_1,D_2,beta=M_1,iter=j,reg_type = 'poisson',formula = f)
  if(M_2$convergence)break
}
Summary_FL_GLM(M_2)

summary(m <- glm(num_awards ~ prog + math, family="poisson", data=p))


#####binomial example

df <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
df$rank <- as.factor(df$rank)

df1=df[1:(nrow(df)/2),]
df2=df[-(1:(nrow(df)/2)),]
maxit=25 
tol=1e-16
M_1=NULL
f=admit ~ gre+gpa+rank
#iterations
for(j in 1:maxit){
  N_1=node_beta(formula =f,Data = df1,beta = M_1$coef,iter = j,reg_type = 'binomial')
  N_2=node_beta(formula = f,Data = df2,beta = M_1$coef,iter = j,reg_type = 'binomial')
  M_1=master_beta(N_1,N_2,beta = M_1,reg_type = 'binomial')
  D_1=node_deviance(formula = f,Data = df1,beta = M_1$coef,iter = j,reg_type = 'binomial')
  D_2=node_deviance(formula = f,Data = df2,beta = M_1$coef,iter = j,reg_type = 'binomial')
  M_2=master_convergence(D_1,D_2,beta=M_1,iter=j,reg_type = 'binomial',formula = f)
  if(M_2$convergence)break
}
Summary_FL_GLM(M_2)
summary(logit <- glm(admit ~ gre+gpa+rank,data=df,family="binomial"))



#####gaussian example

df <- anorexia

df1=anorexia[1:(nrow(anorexia)/2),]
df2=anorexia[-(1:(nrow(anorexia)/2)),]
maxit=25 
tol=1e-16
M_1=NULL
f=Postwt ~ Prewt + Treat + offset(Prewt)
#iterations
for(j in 1:maxit){
  N_1=node_beta(formula =f,Data = df1,beta = M_1$coef,iter = j,reg_type = 'gaussian')
  N_2=node_beta(formula = f,Data = df2,beta = M_1$coef,iter = j,reg_type = 'gaussian')
  M_1=master_beta(N_1,N_2,beta = M_1,reg_type = 'gaussian')
  D_1=node_deviance(formula = f,Data = df1,beta = M_1$coef,iter = j,reg_type = 'gaussian')
  D_2=node_deviance(formula = f,Data = df2,beta = M_1$coef,iter = j,reg_type = 'gaussian')
  M_2=master_convergence(D_1,D_2,beta=M_1,iter=j,reg_type = 'gaussian',formula = f)
  if(M_2$convergence)break
}
Summary_FL_GLM(M_2)
summary(anorex.1 <- glm(Postwt ~ Prewt + Treat + offset(Prewt),family = gaussian, data = anorexia))










#################
p <- read.csv("https://stats.idre.ucla.edu/stat/data/poisson_sim.csv")
p <- within(p, {
  prog <- factor(prog, levels=1:3, labels=c("General", "Academic", 
                                            "Vocational"))
  id <- factor(id)
})
summary(p)

with(p, tapply(num_awards, prog, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))

summary(m <- glm(num_awards ~ prog + math, family="poisson", data=p))


f=summary(m)
f$cov.unscaled

mm=m$fitted.values
fitted(m)
cbind(mu,fitted(m))
sqrt(diag(summary(m)$cov.unscaled)*summary(m)$dispersion)
W=diag(mm*(1-mm)) 
solve(t(X)%*%W%*%X)
vcov(m)

W=diag(as.vector(mu)*(1-as.vector(mu)))
V=solve(t(X)%*%W%*%X)
diag(as.vector(mu*(1-mu)))
sqrt(diag(V))

#######

utils::data(anorexia, package = "MASS")

anorex.1 <- glm(Postwt ~ Prewt + Treat ,
                family = gaussian, data = anorexia)
summary(anorex.1)
X <- model.matrix(Postwt ~ Prewt + Treat ,data = anorexia)
glm_irls(X,anorexia$Postwt,family = gaussian)

