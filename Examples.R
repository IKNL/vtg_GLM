remove(list=ls())
source("Poisson_FL.R")
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
  M_2=master_convergence(D_1,D_2,beta=M_1,iter=j,reg_type = 'poisson',formula = f,maxit=maxit)
  if(M_2$convergence)break
}
Summary_FL_GLM(M_2)

summary(m <- glm(num_awards ~ prog + math, family="poisson", data=p))
#summary(m <- glm(num_awards ~ 1, family="poisson", data=p))


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
  M_2=master_convergence(D_1,D_2,beta=M_1,iter=j,reg_type = 'binomial',formula = f,maxit=maxit)
  if(M_2$convergence)break
}
Summary_FL_GLM(M_2)
summary(logit <- glm(admit ~ gre+gpa+rank,data=df,family="binomial"))



#####gaussian example

df <- anorexia

df1=anorexia[1:(nrow(anorexia)/2),]
df2=anorexia[-(1:(nrow(anorexia)/2)),]
maxit=35 
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
  M_2=master_convergence(D_1,D_2,beta=M_1,iter=j,reg_type = 'gaussian',formula = f,maxit=maxit,tol=tol)
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

