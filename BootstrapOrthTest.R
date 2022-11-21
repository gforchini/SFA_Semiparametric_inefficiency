library(mgcv)
library(stargazer)
library(GA)
library(compiler)
library(fastDummies)
library(doParallel)
library(latex2exp)
library(plotrix)
library(foreach)
library(doRNG)
library(snow)
library(data.table)

######################################################################
############ DGP Simulations Function ################################
######################################################################


# This function simulate the data for the different design 
# defined by lambda and k
Model<-cmpfun(function(n,alpha,beta,sigmaU,sigmaV2,k){
  # DESCRIBE THE INPUTS' DIMENSIONS
  # Simulates the covariates
  X1<-rnorm(n,4,2)
  X2<-rnorm(n,4,2)
  X3<-.5*X2+.5*rnorm(n,4,2)
  X4<-rnorm(n,4,2)
  x<-cbind(X1,X2,X3,X4)
  
  
  ## Simulates the inefficiency term
  U<-abs(rnorm(n,sd=k*sigmaU(x%*%alpha)))
  
  ## Simulates the error term V
  V<-rnorm(n,sd=sqrt(sigmaV2))
  
  ## Simulates the explained variable
  y<-1+x%*%beta-U+V
  
  return(list(y,x)) # Return the simulated data
})


# defines alpha in terms of phi
# the output is a vector
# maps a k-1 dimensional vector into a k dimensional vector on the unit sphere.
A<-cmpfun(function(phi){ c(1,phi)/sqrt(1+sum(phi^2))})

# Jacobian of phi->alpha. Code from Wood
# The output is a matrix of dimension k=length(alpha) by k-1=length(phi)
J<-cmpfun(function(phi){
  kk<-1+sum(phi^2)
  J <- outer(A(phi),-phi/(kk)) ## compute Jacobian
  for (j in 1:length(phi) ){ J[j+1,j] <- J[j+1,j] + 1/sqrt(kk)}
  return(J)
})


###################### estimation accounting for heteroskedasticity ###################3
SP <- cmpfun(function(gamma,y,x,w, opt=TRUE) {
  # This function re-estimate alpha and beta with transformed data that deal
  # with heteroscedasticity
  # y is a vector of dimension n x 1
  # x is a matrix of dimension n x p
  # gamma is a vector of dimension 2k-1 x 1. It contains alpha p x 1 and beta p x 1.
  # sig is a vector of estimated variances for each i
  # Return ML if opt=TRUE and fitted gam with theta added otherwise.
  
  p<-dim(x)[2] # Number of rows of x. This is the dimension of alpha and beta
  phi<-gamma[1:(p-1)]  # initial values for phi
  beta<-gamma[p:(2*p-1)] # initial values for beta
  alpha<-A(phi) # map phi to alpha
  
  a <- x%*%alpha 
  b <- x%*%beta
  y2<-(b-y)*w 
  # w is the weight that corrects for heteroscedasticity
  e<- gam(y2~-1+w+s(a,bs="cr",by=w),method="REML") ## fit the varying-coefficient model
  
  if (opt) 
    return(e$gcv.ubre) 
  else {
    e$phi <- phi
    e$alpha <- alpha
    e$beta <- beta
    e$regressors <- x
    return(e)
  }
})



###################### estimation accounting for heteroskedasticity in bootstrap ###################3
SPb <- cmpfun(function(gamma,y,x,w, sp,opt=TRUE) {
  # This function re-estimate alpha and beta with transformed data that deal
  # with heteroscedasticity
  # y is a vector of dimension n x 1
  # x is a matrix of dimension n x p
  # gamma is a vector of dimension 2k-1 x 1. It contains alpha p x 1 and beta p x 1.
  # sig is a vector of estimated variances for each i
  # Return ML if opt=TRUE and fitted gam with theta added otherwise.
  
  p<-dim(x)[2] # Number of rows of x. This is the dimension of alpha and beta
  phi<-gamma[1:(p-1)]  # initial values for phi
  beta<-gamma[p:(2*p-1)] # initial values for beta
  alpha<-A(phi) # map phi to alpha
  
  a <- x%*%alpha 
  b <- x%*%beta
  y2<-(b-y)*w 
  # w is the weight that corrects for heteroscedasticity
  e<- gam(y2~-1+w+s(a,bs="cr",by=w),sp=sp) ## fit the varying-coefficient model
  
  if (opt) 
    return(e$gcv.ubre) 
  else {
    e$phi <- phi
    e$alpha <- alpha
    e$beta <- beta
    e$regressors <- x
    return(e)
  }
})


#############################################################################################
####### Wild bootstrap without constraints on alpha and heteroscedasticity correction #######
#############################################################################################

boot_stat<-cmpfun(function(fit){
  # input is fitted model
  # output is boostrapped orthotonality test statistic
  eps<-rnorm(n)
  X<-fit$regressors # Extract regressors
  n<-nrow(X)
  # set up Bootstrap model
  beta_h<-fit$beta
  phi_h<-fit$phi
  p_ini<-c(phi_h,beta_h)
  lambda_h<-fit$fitted.values
  s<-fit$residuals
  y_b<-X%*%beta_h-lambda_h+s*eps  # Bootstrap model
  # fit bootstrap model
  fit_b<-optim(p_ini,SP,y=y_b, x=X,w=rep(1,n),method = "BFGS")
  f_b<-SPb(fit_b$par,y=y_b,x=X,w=rep(1,n),fit$sp,opt=FALSE) ## extract best fit model
  alpha_hat<-A(f_b$phi)
  beta_hat<-f_b$beta
  
  return( sqrt(200)*sum(alpha_hat*beta_hat) ) # return alpha'beta
})
  

bootstrap_test_statistic<-cmpfun(function(fit){
  rep_b<-100
  results<-replicate(rep_b,boot_stat(fit))
  s2<-var(results)
  t<-200*(sum((A(fit$phi)*fit$beta))^2)/s2
  return(t)
})

################################################################
##### now computation
################################################################

Estimation<-cmpfun(function(y,x,k,sigmaV2){
  
  # Set starting values close to true value of parameters to speed up convergence
  phi<-alpha[2:length(alpha)]/alpha[1]
  eps_phi<-0.1*phi*runif(3,-1,1)# Deviation between the initial value and the true values
  eps_beta<-0.1*beta*runif(4,-1,1)
  eps_sigmaV2<-0.1*sqrt(sigmaV2)*runif(1,-1,1)
  eps_k<-0.1*sqrt(k)*runif(1,-1,1)
  eps_eta<-0.1*1*runif(1,-1,1)
  
  
  initial_values_sp<-c(phi+eps_phi,beta+eps_beta) # Initial values for SP
  
  # SP Estimation without heteroskedasticity correction
  f1 <-optim(initial_values_sp, SP, x=x,y=y,w=rep(1,n),method = "BFGS")  #"Nelder-Mead"
  fit1 <- SP(f1$par,y=y,x=x,w=rep(1,n),opt=FALSE) # Return the estimates
  
  
  
  test_orth<-bootstrap_test_statistic(fit1)
  
  
  return( test_orth)
  
})




Simulation<-function(alpha, beta,sigmaU,rep,n,titletable){
  # Exact value of Var(X'beta) 
  Var_xbeta<-(beta[1]^2+(beta[2]+.5*beta[3])^2+(.5*beta[3])^2+beta[4]^2)*4
  
  # Numerical approximation of k 
  set.seed(5521249)
  X1=rnorm(10000000,4,2)
  X2=rnorm(10000000,4,2)
  X3=.5*X2+.5*rnorm(10000000,4,2)
  X4=rnorm(10000000,4,2)
  x=cbind(X1,X2,X3,X4)
  
  U<-abs(rnorm(10000000,sd=sigmaU(x%*%alpha)))
  k1<-as.numeric(sqrt((.5*Var_xbeta/1.1)/var(U)))
  
  U<-abs(rnorm(10000000,sd=k1*sigmaU(x%*%alpha)))
  delta1<-mean(U)
  eta<-1
  
  mu1<-eta-delta1
  sigmaV2<-0.1*var(U)
  
  Design<-foreach (i=1:rep,.export=c("gam","Model","Estimation","SP","A","J",
                                     "mu1","bootstrap_test_statistic","boot_stat","SPb"),.combine=cbind,.options.RNG = 58219)  %dorng% {
    
    S <-Model(n=n,alpha=alpha,beta=beta,sigmaU=sigmaU,sigmaV2=sigmaV2,k=k1)
    
    ########### Extract sample ##########
    y<-S[[1]] 
    x<-S[[2]]
    
    
    ########## Run Estimations ###########
    Estimation(y=y,x=x,k=k1,sigmaV2=sigmaV2)
  }
  
  
  t<-Reduce(cbind,Design[1,])
  
  
  Design.10<-mean(t>2.706)
  Design.5<-mean(t>3.841)
  Design.1<-mean(t>6.635)
  
  print(alpha)
  print(beta)
  print(paste("10% rejection probability: ",Design.10))
  print(paste("5% rejection probability: ",Design.5))
  print(paste("1% rejection probability: ",Design.1))
  write.csv(data.frame(t),paste(titletable,"_orth.csv"))
  #return(t)
}





rep<-1000
n<-200
workers <- 12
cl <- parallel::makeCluster(workers)
registerDoParallel(cl)
############################################################
# Design 5: trigonometric with separability 
############################################################

alpha<-c(.6,-.8,0 ,0)
beta<-c(0,0,-.3,.8)
sigmaU<- function(x) 1+cos(x)
titletable<-"Design_5"

Simulation(alpha, beta,sigmaU,rep,n,titletable)


###############################################################
###### Design 6: trigonometric without separability ###########
###############################################################

alpha<-c(0.5, -0.5, -0.5,  0.5)
beta<-c(.75,-.25,-.25,.75)
sigmaU<- function(x) 1+cos(x)
titletable<-"Design 6"

Simulation(alpha, beta,sigmaU,rep,n,titletable)

############################################################
###### Design 7: locally linear with separability ##########
############################################################

alpha<-c(.6,-.8,0 ,0)
alpha<-alpha/as.vector(sqrt(alpha%*%alpha)) #Alpha is on the unit sphere 
beta<-c(0,0,-.3,.8)
sigmaU<- function(x) {ifelse(( (x+.5)/2 <= -1),25+22*(x+.5)+5*(x+.5)^2,ifelse((-1<(x+.5)/2),5+2*(x+.5),NA))}
titletable<-"Design 7"

Simulation(alpha, beta,sigmaU,rep,n,titletable)


############################################################
###### Design 8: locally linear without separability #######
############################################################

alpha<-c(0.5, -0.5, -0.5,  0.5)
alpha<-alpha/as.vector(sqrt(alpha%*%alpha)) #Alpha is on the unit sphere 
beta<-c(.75,-.25,-.25,.75)
sigmaU<- function(x) {ifelse(( (x+.5)/2 <= -1),25+22*(x+.5)+5*(x+.5)^2,ifelse((-1<(x+.5)/2),5+2*(x+.5),NA))}
titletable<-"Design 8"

Simulation(alpha, beta,sigmaU,rep,n,titletable)

############################################################
###### Design 9: partially constant with separability ######
############################################################

alpha<-c(.6,-.8,0 ,0)
alpha<-alpha/as.vector(sqrt(alpha%*%alpha)) #Alpha is on the unit sphere 
beta<-c(0,0,-.3,.8)

sigmaU<- function(x) {ifelse(( x < -4*pi/9 ),1+cos(2.25*x),ifelse(( x > 4*pi/9 ),1+cos(2.25*x), 0))}
titletable<-"Design 9"

Simulation(alpha, beta,sigmaU,rep,n,titletable)


################################################################
###### Design 10: partially constant without separability ######
################################################################

alpha<-c(0.5, -0.5, -0.5,  0.5)
alpha<-alpha/as.vector(sqrt(alpha%*%alpha)) #Alpha is on the unit sphere 
beta<-c(.75,-.25,-.25,.75)
sigmaU<- function(x) {ifelse(( x < -4*pi/9 ),1+cos(2.25*x),ifelse(( x > 4*pi/9 ),1+cos(2.25*x), 0))}
titletable<-"Design 10"

Simulation(alpha, beta,sigmaU,rep,n,titletable)


######################################################
#The design 1 has exponential Inefficiency and separability
######################################################
alpha<-c(.6,-.8,0 ,0)
beta<-c(0,0,-.3,.8)
sigmaU<-exp
titletable<-"Design 1"

Simulation(alpha, beta,sigmaU,rep,n,titletable)


######################################################
# Design 2: Logistic with separability 
######################################################

alpha<-c(.6,-.8,0 ,0)
beta<-c(0,0,-.3,.8)
sigmaU<-function(x) exp(x)/(1+exp(x))
titletable<-"Design 2"

Simulation(alpha, beta,sigmaU,rep,n,titletable)

#######################################
# Design 3: Exponential
#######################################

alpha<-c(0.5, -0.5, -0.5,  0.5)
beta<-c(.75,-.25,-.25,.75)
lambda<-exp
titletable<-"Design 3"

Simulation(alpha, beta,sigmaU,rep,n,titletable)



#################################
# Design 4: Logistic 
#################################


alpha<-c(0.5, -0.5, -0.5,  0.5)
beta<-c(.75,-.25,-.25,.75)
sigmaU<-function(x) exp(x)/(1+exp(x))
titletable<-"Design 4"

Simulation(alpha, beta,sigmaU,rep,n,titletable)


