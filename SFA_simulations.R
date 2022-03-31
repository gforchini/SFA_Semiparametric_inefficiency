
#################################################################
## This files contains the functions used in "Simulations SFA" ##
#################################################################



##################################
### Load the required packages ###
##################################

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


####################################################################################
############ Estimator Functions ################################
####################################################################################


############################## Log-likelihood for ML ############
MLE <- cmpfun(function(gamma,y,x,lambda,opt=TRUE) {
  # DESCRIBE THE INPUTS' DIMENSIONS
  p<-dim(x)[2] # Number of rows of x. This is the dimension of alpha and beta
  alpha<-gamma[1:p] # the vector of the first p components of gamma is alpha
  beta<-gamma[(p+1):(2*p)] # the vector of the last p components of gamma is beta
  sigmaV<-gamma[2*p+1] # Variance of V
  k<-gamma[2*p+2] # Parameter that adapts the parametric model to the
  eta<-gamma[2*p+3] # intercept
  # scale of inefficient
  a <- x%*%alpha 
  b <- x%*%beta 
  delta<-(k^2)*lambda(a)/sigmaV
  r<-y-b- eta
  sigma<-sqrt((k^4)*lambda(a)^2+sigmaV^2)
  #Log-likelihood function
  llik <- sum(-log(sigma) - (1/(2*sigma^2))* (r^2) + pnorm(-r*delta/sigma,0,1, log.p=TRUE))
  if (opt) return(-llik) else {
    return(c(alpha,beta,k,eta))  ## add alpha, beta, k and sigmav^2
  }
})

################# log-likelihood for ML with Separability Constrain ############
MLESep <- cmpfun(function(gamma,y,x,lambda,opt=TRUE) {
  # DESCRIBE THE INPUTS' DIMENSIONS
  p<-dim(x)[2] # Number of rows of x. This is the dimension of alpha and beta
  alpha<-c(gamma[1:(p/2)],rep(0,p/2)) # last two coefficients constrained to 0
  beta<-c(rep(0,p/2),gamma[(p+p/2+1):(2*p)]) # two first coefficients constrained to 0
  sigmaV<-gamma[2*p+1]
  k<-gamma[2*p+2]
  eta<-gamma[2*p+3] # intercept
  a <- x%*%alpha 
  b <- x%*%beta 
  delta<-(k^2)*lambda(a)/sigmaV
  r<-y-b-eta #Remove the linear part
  sigma<-sqrt((k^4)*lambda(a)^2+sigmaV^2)
  llik <- sum(-log(sigma) - (1/(2*sigma^2))* (r^2) + pnorm(-r*delta/sigma,0,1, log.p=TRUE))
  if (opt) return(-llik) else {
    return(c(alpha,beta,k,eta))  
  }
})




###################### estimation accounting for heteroskedasticity ###################3
SP <- cmpfun(function(gamma,y,x,w=w, opt=TRUE) {
  # This function re-estimate alpha and beta with tranformed data that deal
  # with heteroscedasticity
  # y is a vector of dimension n x 1
  # x is a matrix of dimension n x p
  # gamma is a vector of dimension 2k x 1. It contains alpha p x 1 and beta p x 1.
  # sig is a vector of estimated variances for each i
  # Return ML if opt=TRUE and fitted gam with theta added otherwise.
  
  p<-dim(x)[2] # Number of rows of x. This is the dimension of alpha and beta
  alpha<-gamma[1:p] # the vector of the first p components of gamma is alpha
  alpha[1]<-abs(alpha[1]) # The first component of alpha is positive (identification)
  beta<-gamma[(p+1):(2*p)] 
  kk <- sum(alpha^2) # sum of squares of alpha
  alpha <- alpha/sqrt(kk)  ## normalized alpha 
  a <- x%*%alpha 
  b <- x%*%beta 
  y2<-(b-y)*w 
  # w is the weight that corrects for heteroscedasticity
  e<- gam(y2~-1+w+s(a,bs="cr",by=w),method="REML") ## fit the varying-coefficient model
  
  if (opt) return(e$gcv.ubre) else {
    e$alpha <- alpha 
    e$beta <- beta           
    return(e)
  }
})

####################################################################################
############ Estimation Function ################################
####################################################################################



Estimation<-cmpfun(function(y,x,k,sigmaV2){
  
  # EXPLAIN WHAT IS HAPPENING IN THE LINES BELOW
  
  eps_alpha<-0.1*alpha*runif(4,-1,1)# Deviation between the initial value and the true values
  eps_beta<-0.1*beta*runif(4,-1,1)
  eps_sigmaV2<-0.1*sqrt(sigmaV2)*runif(1,-1,1)
  eps_k<-0.1*sqrt(k)*runif(1,-1,1)
  eps_eta<-0.1*1*runif(1,-1,1)
  
  
  initial_values<-c(alpha+eps_alpha,beta+eps_beta,sqrt(sigmaV2)+eps_sigmaV2,sqrt(k)+eps_k,1+eps_eta) # Initial values
  
  
  # Estimation without heteroscedasticity correction
  f1 <-optim(initial_values[1:8], SP, x=x,y=y,w=rep(1,n),method = "BFGS") 
  fit1 <- SP(f1$par,y=y,x=x,w=rep(1,n),opt=FALSE)# Return the estimates
  
  # ExtraCtion of the heterosedastic residuals
  a<-x%*%fit1$alpha
  b<-x%*%fit1$beta
  
  y2<-b-y  
  ## fit model
  v2<-log(fit1$residuals^2)
  e<-gam(v2~s(a,bs="cr"),method="REML")
  
  
  w<-1/sqrt(exp(e$fitted.values)) #weights that correct for heteroscedasticity
  
  f2 <-optim(c(fit1$alpha,fit1$beta),SP,x=x,y=y,w=w,method = "BFGS")
  
  
  #Estimation with heteroscedasticiy correction
  fit2 <- SP(f2$par,x=x,y=y,w=w,opt=FALSE) ## extract best fit model

  
  
  a<-x%*%fit2$alpha
  b<-x%*%fit2$beta
  
  
  ##### General MLE ######
  
  # MLE with exponential inefficiencies
  f0MLE <-optim(initial_values,MLE,lambda=exp, x=x,y=y,method = "BFGS") 
  theta.est<-f0MLE$par
  fitMLEexp <- MLE(theta.est,lambda=exp, x=x,y=y,opt=FALSE)
  
  
  
  # MLE with logistic inefficiencies
  f0MLE <-optim(initial_values,MLE,lambda=plogis, x=x,y=y,method = "BFGS") 
  theta.est<-f0MLE$par
  fitMLElogis <- MLE(theta.est,lambda=plogis, x=x,y=y,opt=FALSE)
  
  
  ##### MLE with separability #####
  
  # Note that the Nelder-Mead method is used  instead of BFGS for the exponential
  
  # MLE with exponential inefficiencies
  f0MLE <-optim(initial_values,MLESep,lambda=exp, x=x,y=y,method = "Nelder-Mead") 
  
  
  theta.est<-f0MLE$par
  fitMLESepexp <- MLESep(theta.est,lambda=exp, x=x,y=y,opt=FALSE)
  
  
  # MLE with logistic inefficiencies
  f0MLE <-optim(initial_values,MLESep,lambda=plogis, x=x,y=y,method = "Nelder-Mead") 
  
  
  theta.est<-f0MLE$par
  fitMLESeplogis <- MLESep(theta.est,lambda=plogis, x=x,y=y,opt=FALSE)
  
  
  # Estimation of the intercepts
  delta.MLEexp<-sqrt(2/pi)*mean(fitMLEexp[9]^2*exp(x%*%fitMLEexp[1:4]))
  delta.MLElogis<-sqrt(2/pi)*mean(fitMLElogis[9]^2*plogis(x%*%fitMLElogis[1:4]))
  delta.MLESepexp<-sqrt(2/pi)*mean(fitMLESepexp[9]^2*exp(x%*%fitMLESepexp[1:4]))
  delta.MLESeplogis<-sqrt(2/pi)*mean(fitMLESeplogis[9]^2*plogis(x%*%fitMLESeplogis[1:4]))
  delta.hom<--min(predict(fit1,type = "terms")[,2])
  delta.hete<-- min(predict(fit2,type="terms")[,2]/w)
  
  
  
  eta.MLEexp<-fitMLEexp[10]
  eta.MLElogis<-fitMLElogis[10]
  eta.MLESepexp<-fitMLESepexp[10]
  eta.MLESeplogis<-fitMLESeplogis[10]
  
  mu.MLEexp<-eta.MLEexp-delta.MLEexp
  mu.MLElogis<-eta.MLElogis-delta.MLElogis
  mu.MLESepexp<-eta.MLEexp-delta.MLESepexp
  mu.MLESeplogis<-eta.MLEexp-delta.MLESeplogis
  mu.hom<-mean(y-x%*%fit1$beta)
  mu.hete<-mean(y-x%*%fit2$beta)
  
  
  
  
  eta.hom<- mu.hom+delta.hom
  eta.hete<-mu.hete+delta.hete
  
  # Return the estimates
  results<-list(c(delta.MLEexp,
                  delta.MLElogis,
                  delta.MLESepexp,
                  delta.MLESeplogis,
                  delta.hom,
                  delta.hete),
                c(mu.MLEexp,
                  mu.MLElogis,
                  mu.MLESepexp,
                  mu.MLESeplogis,
                  mu.hom,
                  mu.hete),
                c(eta.MLEexp,
                  eta.MLElogis,
                  eta.MLESepexp,
                  eta.MLESeplogis,
                  eta.hom,
                  eta.hete),
                c(fitMLEexp[1:4],fitMLElogis[1:4],
                  fitMLESepexp[1:4],fitMLESeplogis[1:4],
                  fit1$alpha,fit2$alpha),
                c(fitMLEexp[5:8],fitMLElogis[5:8],
                  fitMLESepexp[5:8],fitMLESeplogis[5:8],
                  fit1$beta,fit2$beta))
  
  return( results)
  
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
  
  Design<-foreach (i=1:rep,.export=c("gam","Model","Estimation","SP","MLE","MLESep","mu1"),.combine=cbind,.options.RNG = 58219)  %dorng% {
    
    S <-head(Model(n=n,alpha=alpha,beta=beta,sigmaU=sigmaU,sigmaV2=sigmaV2,k=k1))
   
    ########### Extract sample ##########
    y<-S[[1]] 
    x<-S[[2]]
    
    
    ########## Run Estimations ###########
    Estimation(y=y,x=x,k=k1,sigmaV2=sigmaV2)
    }
  
  
  deltaDesign<-Reduce(cbind,Design[1,])
  muDesign<-Reduce(cbind,Design[2,])
  etaDesign<-Reduce(cbind,Design[3,])
  alphaDesign<-Reduce(cbind,Design[4,])
  betaDesign<-Reduce(cbind,Design[5,])
  

  
  
  
  # Median of the delta, mu and eta
  Design.deltamedian<-apply(deltaDesign,1,median)
  Design.mumedian<-apply(muDesign,1,median)
  Design.etamedian<-apply(etaDesign,1,median)
  
  
  
  # Median of alpha's
  Design.alpha1median<-apply(alphaDesign,1,median)[c(1,5,9,13,17,21)]
  Design.alpha2median<-apply(alphaDesign,1,median)[1+c(1,5,9,13,17,21)]
  Design.alpha3median<-apply(alphaDesign,1,median)[2+c(1,5,9,13,17,21)]
  Design.alpha4median<-apply(alphaDesign,1,median)[3+c(1,5,9,13,17,21)]
  
  
  # Median For beta's
  Design.beta1median<-apply(betaDesign,1,median)[c(1,5,9,13,17,21)]
  Design.beta2median<-apply(betaDesign,1,median)[1+c(1,5,9,13,17,21)]
  Design.beta3median<-apply(betaDesign,1,median)[2+c(1,5,9,13,17,21)]
  Design.beta4median<-apply(betaDesign,1,median)[3+c(1,5,9,13,17,21)]
  
  # Median biases
  Bias.Median.delta.Design<-Design.deltamedian-delta1
  Bias.Median.eta.Design<-Design.etamedian-eta
  Bias.Median.mu.Design<-Design.mumedian-mu1
  Bias.Median.Alpha.Design<-cbind(Design.alpha1median-alpha[1],Design.alpha2median-alpha[2],Design.alpha3median-alpha[3],Design.alpha4median-alpha[4])
  Bias.Median.Beta.Design<-cbind(Design.beta1median-beta[1],Design.beta2median-beta[2],Design.beta3median-beta[3],Design.beta4median-beta[4])
  
  
  
  ##### Compute the MAD of the estimates for each sample #####
  MAD.delta.Design<-apply(abs(deltaDesign-Design.deltamedian),1,median)
  MAD.mu.Design<-apply(abs(muDesign-Design.mumedian),1,median)
  MAD.eta.Design<-apply(abs(etaDesign-Design.etamedian),1,median)
  
  
  MAD.Alpha.Design<-cbind(apply(abs(alphaDesign[c(1,5,9,13,17,21),]-Design.alpha1median),1,median)
                          ,apply(abs(alphaDesign[c(1,5,9,13,17,21)+1,]-Design.alpha2median),1,median)
                          ,apply(abs(alphaDesign[c(1,5,9,13,17,21)+2,]-Design.alpha3median),1,median)
                          ,apply(abs(alphaDesign[c(1,5,9,13,17,21)+3,]-Design.alpha4median),1,median))
  
  MAD.Beta.Design<-cbind(apply(abs(betaDesign[c(1,5,9,13,17,21),]-Design.beta1median),1,median)
                         ,apply(abs(betaDesign[c(1,5,9,13,17,21)+1,]-Design.beta2median),1,median)
                         ,apply(abs(betaDesign[c(1,5,9,13,17,21)+2,]-Design.beta3median),1,median)
                         ,apply(abs(betaDesign[c(1,5,9,13,17,21)+3,]-Design.beta4median),1,median))
  
  ##### Build tables of results #####
  
  # For different intercepts
  
  True.Design.intercept<-t(c(delta1,eta,mu1,0,0,0))
  Result.Design.intercept<-rbind(True.Design.intercept,cbind(Bias.Median.delta.Design,Bias.Median.eta.Design,Bias.Median.mu.Design,MAD.delta.Design,MAD.mu.Design,MAD.eta.Design))
  rownames(Result.Design.intercept)<-c("True","ML(e)","ML(l)","ML_s(e)","ML_s(l) ","PS","WPS")
  colnames(Result.Design.intercept)<-c("MB.delta","MB.eta","MB.mu",
                                       "MAD.delta","MAD.eta","MAD.mu")
  
  
  # For alpha's
  
  Result.Design.alpha<-data.frame(cbind(Bias.Median.Alpha.Design,MAD.Alpha.Design))               
  rownames(Result.Design.alpha)<-c("ML(e)","ML(l)","ML_s(e)","ML_s(l) ","PS","WPS")
  colnames(Result.Design.alpha)<-c("MB.alpha1","MB.alpha2","MB.alpha3","MB.alpha4",
                                   "MAD.alpha1","MAD.alpha2","MAD.alpha3","MAD.alpha4")
  
  
  # For beta's
  
  Result.Design.Beta<-data.frame(cbind(Bias.Median.Beta.Design,MAD.Beta.Design))
  rownames(Result.Design.Beta)<-c("ML(e)","ML(l)","ML_s(e)","ML_s(l) ","PS","WPS")
  colnames(Result.Design.Beta)<-c("MB.beta1","MB.beta2","MB.beta3","MB.beta4",
                                  "MAD.beta1","MAD.beta2","MAD.beta3","MAD.beta4")
  
  
  # Tex format
  stargazer(Result.Design.intercept,digits = 3,summary = FALSE,title=titletable)
  stargazer(Result.Design.alpha,digits = 3,summary = FALSE,title=titletable)
  stargazer(Result.Design.Beta,digits = 3,summary = FALSE,title=titletable)
  
}




#########################################################################
# Set sample size, number of replication and initialize cluster for parallel
# computations
#########################################################################
n<-200  # Sample size
rep<-10000 # Number of samples

workers <- 27
cl <- parallel::makeCluster(workers)

# register the cluster for using foreach
doParallel::registerDoParallel(cl)



######################################################
#The design 1 has exponential Inefficiency and separability
######################################################
alpha<-c(.6,-.8,0 ,0)
alpha<-alpha/as.vector(sqrt(alpha%*%alpha)) #Alpha is on the unit sphere 
beta<-c(0,0,-.3,.8)
sigmaU<-exp
titletable<-"Design 1"

Simulation(alpha, beta,sigmaU,rep,n,titletable)


######################################################
# Design 2: Logistic with separability 
######################################################

alpha<-c(.6,-.8,0 ,0)
alpha<-alpha/as.vector(sqrt(alpha%*%alpha)) #Alpha is on the unit sphere 
beta<-c(0,0,-.3,.8)
sigmaU<-function(x) exp(x)/(1+exp(x))
titletable<-"Design 2"

Simulation(alpha, beta,sigmaU,rep,n,titletable)

#######################################
# Design 3: Exponential
#######################################

alpha<-c(1,-1, -1 , 1)
alpha<-alpha/as.vector(sqrt(alpha%*%alpha)) #Alpha is on the unit sphere 
beta<-c(.75,-.25,-.25,.75)
lambda<-exp
titletable<-"Design 3"

Simulation(alpha, beta,sigmaU,rep,n,titletable)



#################################
# Design 4: Logistic 
#################################


alpha<-c(1,-1, -1 , 1)
alpha<-alpha/as.vector(sqrt(alpha%*%alpha)) #Alpha is on the unit sphere 
beta<-c(.75,-.25,-.25,.75)
sigmaU<-function(x) exp(x)/(1+exp(x))
titletable<-"Design 4"

Simulation(alpha, beta,sigmaU,rep,n,titletable)



############################################################
# Design 5: trigonometric with separability 
############################################################

alpha<-c(.6,-.8,0 ,0)
alpha<-alpha/as.vector(sqrt(alpha%*%alpha)) #Alpha is on the unit sphere 
beta<-c(0,0,-.3,.8)
sigmaU<- function(x) 1+cos(x)
titletable<-"Design 5"

Simulation(alpha, beta,sigmaU,rep,n,titletable)


###############################################################
###### Design 6: trigonometric without separability ###########
###############################################################

alpha<-c(1,-1, -1 , 1)
alpha<-alpha/as.vector(sqrt(alpha%*%alpha)) #Alpha is on the unit sphere 
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

alpha<-c(1,-1, -1 , 1)
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

alpha<-c(1,-1, -1 , 1)
alpha<-alpha/as.vector(sqrt(alpha%*%alpha)) #Alpha is on the unit sphere 
beta<-c(.75,-.25,-.25,.75)
sigmaU<- function(x) {ifelse(( x < -4*pi/9 ),1+cos(2.25*x),ifelse(( x > 4*pi/9 ),1+cos(2.25*x), 0))}
titletable<-"Design 10"

Simulation(alpha, beta,sigmaU,rep,n,titletable)


