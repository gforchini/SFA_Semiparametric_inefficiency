###########################################
####### Used Packages #####################
###########################################

library(xtable)
library(ggplot2)
library(mgcv)
library(stargazer)
library(latex2exp)
library(GA)
library(compiler)
library(fastDummies)
library(foreach)
library(doRNG)
library(doParallel)
library(data.table)

###########################################
####### Estimation functions ##############
###########################################

# loglikelihood with exponential inefficiencies
logLik.exp= function(gamma,y,x,opt=TRUE){
  alpha=c(gamma[1],gamma[2],-gamma[1],gamma[3],0,0,0,rep(0,16)) # alpha has 23 components
  beta=gamma[4:26] #  beta has 23 components
  sigmaU=exp(-gamma[27]-x%*%alpha) #standard deviation function
  sigmaV=exp(gamma[28]) # this is a scalar
  
  # sigma and lambda parameters
  sigma=sqrt(sigmaU^2+sigmaV^2)
  lambda= sigmaU/sigmaV
  
  # residuals
  r= y-gamma[29]-x%*%beta;
  
  # change sign of residuals
  r=-r
  
  # log-likelihood: eq (2) of Green 1982, Journal of Econometrics
  llik = -log(1/2)-log(sigma) +dnorm(r/sigma,0,1,log=TRUE) +  pnorm(-r*lambda/sigma,0,1, log.p=TRUE); 
  
  #The loglikelihood is
  c=list(lik= sum(llik))
  
  if (opt) return(c$lik) else {
    c$alpha =alpha  ## add alpha
    c$beta =beta ## add beta
    c$beta0=gamma[29] # constant in beta
    c$alpha0= -gamma[27] # constant in alpha 
    return(c)
  }
}

# loglikelihood with logistic inefficiencies
logLik.logis= function(gamma,y,x,opt=TRUE){
  
  alpha=c(gamma[1],gamma[2],-gamma[1],gamma[3],0,0,0,rep(0,16)) # alpha has 23 components
  beta=gamma[4:26] #  beta has 23 components
  
  sigmaU = 1/(1+exp(gamma[27]+ x%*%alpha)) # This is a vectors 
  sigmaV=exp(gamma[28]) # this is a scalar
  
  # sigma and lambda parameters
  sigma=(sigmaU^2+sigmaV^2)^(1/2)
  
  lambda=(sigmaU/sigmaV);
  
  # residuals
  r= y-gamma[29]- x%*%beta;
  
  # change sign of residuals
  r = -r
  
  # log-likelihood: eq (2) of Green 1982, Journal of Econometrics
  llik = -log(1/2)-log(sigma) +dnorm(r/sigma,0,1,log=TRUE) +  pnorm(-r*lambda/sigma,0,1, log.p=TRUE);
  sum(llik)
  
  #The  loglikelihood is
  c=list(lik= sum(llik))
  
  if (opt) return(c$lik) else {
    c$alpha =alpha  ## add alpha
    c$beta =beta ## add beta
    c$beta0=gamma[29]
    c$alpha0= -gamma[27] # constant in alpha
    return(c)
  }
}

###################### estimation accounting for heteroskedasticity ###################3
hetero.est <- cmpfun(function(upsilon,e,x, opt=TRUE) {
  # This function estimates the heteroscedasticity
  # e is a vector of dimension n x 1
  # x is a matrix of dimension n x k
  # Return ML if opt==TRUE and fitted gam with theta added otherwise.
  x1=x[,1:7] #continous functions are in the index part
  x2=x[,8:23] #Time dummies are in the linear part
  k<-dim(x1)[2] 
  upsilon[1]<-abs(upsilon[1]) # The first component of upsilon is positive (identification)
  kk <- sum(upsilon^2) # sum of squares of upsilon
  upsilon <- upsilon/sqrt(kk)  ## normalized upsilon
  a <- x1%*%upsilon 
  u <- gam(e~x2+s(a,bs="cr"),family=gaussian(link = "log"),method="REML") ## fit the model
  if (opt) return(u$gcv.ubre) else {
    u$upsilon=upsilon
    return(u)
  }
})

# Penalized Cubic Splines Estimation with constrained alpha
SP=cmpfun(function(gamma,y,x,w=rep(1,816),k=10,sign=1,opt=TRUE,fx=FALSE) {
  # Fit single index model using gam call given alpha. 
  # Return ML if opt==TRUE and fitted gam with theta added otherwise.
  # Suitable for calling from 'optim' to find optimal alpha.
  
  alpha=c(-abs(gamma[1]),gamma[2],abs(gamma[1]),gamma[3],0,0,0,rep(0,16)) # The first component of alpha is positive
  # The first and third component have opposite sign
  # Not all variables are included so some of
  # the components of alpha are zero
  
  beta=gamma[4:26] # beta has 23 components
  
  kk = sqrt(sum(alpha^2))
  alpha = alpha/kk         # so now ||alpha||=1
  
  a = x%*%alpha
  b =x%*%beta 
  
  # residuals
  r= (y-b)*w ;
  
  # This removes the linear part.
  
  c = gam(r~-1+w+s(a,bs="cr",by=w,k=k),method="REML") ## fit single index model with REML criterions
  
  if (opt) return(-sign*c$gcv.ubre) else {
    c$alpha = alpha  ## add alpha 
    c$beta = beta ## add beta Again 
    return(c)
  }
})

# Penalized Cubic Splines Estimation with unconstrained alpha
SP_full= cmpfun(function(gamma,y,x,w=rep(1,816),k=10,sign=1,opt=TRUE,fx=FALSE) {
  # Fit single index model using gam call given alpha. 
  # Return ML if opt==TRUE and fitted gam with theta added otherwise.
  # Suitable for calling from 'optim' to find optimal alpha.
  
  alpha=c(-abs(gamma[1]),gamma[2],gamma[3],gamma[4],gamma[5],gamma[6],gamma[7],rep(0,16)) # The first component of alpha is positive
  # Not all variables are included so some of
  # the components of alpha are zero
  beta=gamma[8:30] # beta has 23 components
  
  kk = sqrt(sum(alpha^2))
  alpha = alpha/kk         # so now ||alpha||=1
  
  a = x%*%alpha
  b =x%*%beta 
  
  # residuals
  r= (y-b)*w # This removes the linear part.
  
  c = gam(r~-1+w+s(a,bs="cr",by=w,k=k),method="REML") # fit single index model with REML criterions
  
  if (opt) return(-sign*c$gcv.ubre) else {
    c$alpha = alpha  ## add alpha 
    c$beta = beta ## add beta Again 
    return(c)
  }
})

#########################################
######## Prepare the data ###############
#########################################

# Load the data set from the address
address="data.txt"
data = read.csv(address, sep="")

# name the variables
colnames(data)=c("Ind", "Year","Q" ,  "Y" ,  "P"  , "POP", "AHS","HDD" ,"CDD",  "SDH")
# Logarithm of the input and output variables
ldata=data.frame(log(data[,3:9]))
# Name the logaritm of the variables
colnames(ldata)=c("ln(Q)" ,  "ln(Y)" ,  "ln(P)"  , "ln(POP)", "ln(AHS)","ln(HDD)" ,"ln(CDD)")
# Merge both data
data=cbind.data.frame(data,ldata)
# Create the time dummy variables
data=fastDummies::dummy_cols(data, select_columns = "Year") # add dummies for year

# Name the input variables
varX=c("ln(Y)","ln(P)","ln(POP)","ln(AHS)","ln(HDD)","ln(CDD)","SDH")
# Name the time dummy variables
Ydummy=c("Year_1996","Year_1997","Year_1998",
         "Year_1999","Year_2000", "Year_2001","Year_2002","Year_2003",
         "Year_2004","Year_2005","Year_2006","Year_2007","Year_2008","Year_2009","Year_2010","Year_2011")

# Define the output variable and input variables
y=as.vector(data[,"ln(Q)"]); 
x=as.matrix(data[,varX]);

# De-mean the variables
y=y-mean(y) 
x=sweep(x, 2, apply(x, 2, mean))

# Add time dummies to x
x=cbind(x,as.matrix(data[,Ydummy])) 

###################################################
#### MLE with exponential inefficiencies ###########
###################################################



# Suggestions close to the optimum
suggestion_exp=c(-6.148206000000,1.445908000000,-12.187251000000, 0.236089000000, -0.113605000000, 
                            0.798683000000, -1.469075000000, 0.347588000000, 0.079834000000, 0.005140000000, 
                            0.048606000000, 0.007129000000, -0.033949000000, -0.031903000000, -0.021668000000, 
                            -0.038934000000, -0.064697000000, -0.052713000000, -0.049432000000,-0.048090000000, 
                            -0.089530000000, -0.077848000000, -0.065927000000, -0.072200000000,-0.113551000000, 
                            -0.098241000000, 4.123592000000,-2.554543000000,0.01911112)

ML_exp = ga(type = "real-valued", y=y, x=x, fitness = logLik.exp, lower = suggestion_exp-.1, 
            upper=suggestion_exp+.1,popSize = 200,seed=97654,run=100,
            maxiter=10000000,parallel = TRUE,suggestions =suggestion_exp )




Result_ML_exp=logLik.exp(suggestion_exp,y,x,opt=FALSE) # Extract best fit model     
Result_ML_exp
#ML=874.9509
Result_ML_exp=logLik.exp(summary(ML_exp)$solution,y,x,opt=FALSE) # Extract best fit model     
Result_ML_exp
#ML=874.9509

Result_ML_exp$alpha0
Result_ML_exp$alpha
Result_ML_exp$beta 


################################################
#### MLE with logistic inefficiencies ##########
################################################


# Suggestions close to the optimum
suggestion_logis=c(-7.014168431000,1.576703399000, -14.282466000000,0.237966675000 ,-0.116957888000 ,
                    0.796928590000 ,-1.479550193000 ,0.347236984000 ,0.080293566000 ,0.005155203000 ,
                    0.048680901000 ,0.007165027000 ,-0.034144138000 ,-0.032125123000 ,-0.021738519000 ,
                    -0.038989271000 ,-0.064807371000 ,-0.052726306000 ,-0.049399308000 ,-0.048128831000 ,
                     -0.089512124000 ,-0.077821669000 ,-0.065646345000 ,-0.071979037000 ,-0.113532550000 ,
                     -0.098062029000 ,4.280481519000 ,-2.554140741000,0.02036251)

ML_logis = ga(type = "real-valued", y=y, x=x, fitness = logLik.logis, lower = suggestion_logis-.001, 
              upper=suggestion_logis+.001,popSize = 200,run=100,seed=97654,
              maxiter=10000,parallel = TRUE,suggestions = suggestion_logis)

Result_ML_logis=logLik.logis(suggestion_logis,y,x,opt=FALSE) # Extract best fit model     
Result_ML_logis
#ML=875.9189
Result_ML_logis=logLik.logis(summary(ML_logis)$solution,y,x,opt=FALSE) # Extract best fit model     
Result_ML_logis
#ML= 875.9189



Result_ML_logis$alpha0
Result_ML_logis$alpha
Result_ML_logis$beta

####################################################################
####### Semi-parametric Estimation with constraints on alpha #######
####################################################################

# Suggestions close to the optimum
suggestion_SP=c(-0.4409693, 0.04101152, -1.516389, 0.2453335, -0.1860581, 0.7824441, -1.574586, 0.318795,
                0.07866492, 0.005507147, 0.05226997, 0.01250153, -0.03482731, -0.03049185, -0.01116791,
                -0.0276573, -0.05319287, -0.03504192, -0.02924658, -0.02617147, -0.06511178, -0.04946462,
                -0.02981384, -0.04029411, -0.08164109, -0.06703873)

# Genetic algorithm that searches for a global optimum between lower and upper with suggested values
ptm = proc.time()
opt = ga(type = "real-valued", y=y, x=x, fitness = SP,k=10,pmutation = 0.9, lower = c(-1,-1,-1,rep(-2,23)), 
         upper =c(0,1,1,rep(2,23)),parallel = TRUE,popSize = 200,run=100,seed=97654, maxiter=100000,suggestions=suggestion_SP)
proc.time() - ptm

# Local search of the optimum with the quasi-Newton method
opt2=optim(summary(opt)$solution,SP,sign=-1,y=y, x=x,k=10,method = "BFGS")
# ML=927.7841

semipara = SP(opt2$par,y,x,opt=FALSE) # Extract best fit model

lambda=predict(semipara ,type="terms")[,2] # Predicted values of lambda

###################################################################################################
####### Semi-parametric Estimation with constraints on alpha and heteroscedastic correction #######
###################################################################################################


# Estimate the heteroscedasticity
v2<-abs(semipara$residuals)
hetero <-optim(semipara$alpha[1:7], hetero.est, x=x,e=v2)  
hetero.fit <- hetero.est(hetero$par,x,e=v2,opt=FALSE) # Extract best fit model
#ML=1395.373

w<-1/exp(hetero.fit$fitted.values) # weights that correct for heteroscedasticity

# Search of the optimum with the quasi-Newton method at the neighborhood of the
# values of the estimates without heteroscedasticity
opt2_hete=optim(summary(opt)$solution,SP,sign=-1,y=y, w=w, x=x,k=10,method = "BFGS")
opt2_hete

#ML=981.6815
semipara_hete = SP(opt2_hete$par,y=y,x=x,w=w,opt=FALSE) # Extract best fit model
semipara_hete 

lambda_hete=predict(semipara_hete ,type="terms")[,2]/w # Predicted values of lambda


####################################################################
####### Semi-parametric Estimation without constraints on alpha ####
####################################################################

# Suggestions close to the optimum
suggestion_full=c(-0.371903109, -0.028991713,  0.292276012, -0.639606972, -0.453739803,  0.042969334,
                  0.025035474,  0.196518852, -0.139507524,  0.791292701, -1.440114889,  0.176623890,
                  0.111788716,  0.019900150,  0.048622118,  0.006661658, -0.035903452, -0.033732869,
                  -0.016957978, -0.037429644, -0.060964750, -0.045560400, -0.044684357, -0.048649355,
                  -0.091417432, -0.084948584, -0.071725595, -0.073539734, -0.106700735, -0.103274060)

# Genetic algorithm that searches for a global optimum between lower and upper with suggested values
ptm = proc.time()
optim_full = ga(type = "real-valued", y=y, x=x, fitness = SP_full,k=10,pmutation = 0.9, lower = c(-1,-1,-1,-1,-1,-1,-1,rep(-2,23)), 
                upper =c(0,1,1,1,1,1,1,rep(2,23)),parallel = TRUE,popSize = 500,run=100
                ,seed=97654, maxiter=100000,suggestions=suggestion_full)
proc.time() - ptm

# Local search of the optimum with the quasi-Newton method
optim_full2=optim(summary(optim_full)$solution,SP_full,sign=-1,y=y, x=x,k=10,method = "BFGS")
# ML=1029.936 

semipara_full=SP_full(optim_full2$par,y,x,opt=FALSE) # Extract best fit model

lambda_full=predict(semipara_full ,type="terms")[,2] # Predicted values of lambda


######################################################################################################
####### Semi-parametric Estimation without constraints on alpha and heteroscedastic correction #######
######################################################################################################

# Estimate the heteroscedasticity
v2<-abs(semipara_full$residuals)
hetero_full <-optim(semipara_full$alpha[1:7], hetero.est, x=x,e=v2) 
hetero.fit_full <- hetero.est(hetero_full$par,x,e=v2,opt=FALSE) # Extract best fit model
#ML=1377.969

w_full<-1/exp(hetero.fit_full$fitted.values) # weights that correct for heteroscedasticity

# Local search of the optimum with the quasi-Newton method
optim_full2_hete=optim(summary(optim_full)$solution,SP_full,sign=-1,y=y, x=x,w=w_full,k=10,method = "BFGS")
#ML=1092.105

semipara_full_hete = SP_full(optim_full2_hete$par,y=y,x=x,w=w,opt=FALSE) # Extract best fit model

lambda_full_hete=predict(semipara_full_hete ,type="terms")[,2]/w # Predicted values of lambda


#######################################################
####### Wild bootstrap with constraints on alpha ######
#######################################################

workers <- 10 # Number of splited cores
cl <- parallel::makeCluster(workers)

# register the cluster for using foreach
doParallel::registerDoParallel(cl)

ptm = proc.time()
rep=2 # number of replicates
Bootstrap=foreach (i=1:rep,.export=c("gam"),.combine=cbind,.options.RNG = 58219)  %dorng% {
  
  eps=rnorm(816) # Resample the error term with a normal distribution
  
  # Resample the output
  y_boot=x%*%semipara$beta+semipara$fitted.values+semipara$residuals*eps
  
  # Estimate the model
  optim_boot=optim(opt2$par,SP,sign=-1,y=y_boot, x=x,k=10,method = "BFGS")
  semipara_boot = SP(optim_boot$par,y_boot,x,k=10,opt=FALSE) # Extract best fit model
  lambda_boot=predict(semipara_boot,type="terms")[,2] # Predicted values of lambda
  
  # Save the results
  results=list(c(semipara_boot$alpha,semipara_boot$beta),lambda_boot,mean(-x%*%semipara_boot$beta))
  
}
proc.time() - ptm

stopCluster(cl)

##########################################################################################
####### Wild bootstrap with constraints on alpha and heteroscedasticity correction #######
##########################################################################################

###### Make the core works in parallel ########

workers <- 10 # Number of splited cores
cl <- parallel::makeCluster(workers)

# register the cluster for using foreach
doParallel::registerDoParallel(cl)

ptm = proc.time()
rep=2 # number of replicates
Bootstrap_hete=foreach (i=1:rep,.export=c("gam"),.combine=cbind,.options.RNG = 58219)  %dorng% {
  
  eps=rnorm(816) # Resample the error term with a normal distribution
  
  # Resample the output
  y_boot=x%*%semipara_hete$beta+semipara_hete$fitted.values+semipara_hete$residuals*eps

  # Estimate the model
  optim_boot=optim(opt2$par,SP,sign=-1,y=y_boot, x=x,k=10,method = "BFGS")
  semipara_boot = SP(optim_boot$par,y_boot,x,k=10,opt=FALSE) # Extract best fit model
  lambda_boot=predict(semipara_boot,type="terms") # Predicted values of lambda
  
  # Estimate the heteroscedasticity
  v2<-abs(semipara_boot$residuals)
  hetero_boot <-optim(semipara_boot$alpha[1:7], hetero.est, x=x,e=v2) 
  hetero.fit <- hetero.est(hetero_boot$par,x,e=v2,opt=FALSE) # Extract best fit model
  w_boot<-1/exp(hetero.fit$fitted.values) # Weights that correct for heteroscedasticity
  
  # Estimate the model
  optim_boot_hete=optim(opt2_hete$par,SP,sign=-1,y=y_boot, x=x,w=w_boot,k=10,method = "BFGS")
  semipara_boot_hete = SP(optim_boot_hete$par,y=y_boot,x=x,w=w_boot,k=10,opt=FALSE) # Extract best fit model
  lambda_boot_hete=predict(semipara_boot_hete,type="terms")[,2]/w_boot # Predicted values of lambda

  # Save the results
  results=list(c(semipara_boot_hete$alpha,semipara_boot_hete$beta),lambda_boot_hete,mean(-x%*%semipara_boot_hete$beta))
  
}
proc.time() - ptm

stopCluster(cl)



##########################################################
####### Wild bootstrap without constraints on alpha ######
##########################################################

###### Make the core works in parallel ########
 
workers <- 10 # Number of splited cores
cl <- parallel::makeCluster(workers)

# register the cluster for using foreach
doParallel::registerDoParallel(cl) 

ptm = proc.time()
rep=2 # number of replicates
Bootstrap_full=foreach (i=1:rep,.export=c("gam"),.combine=cbind,.options.RNG = 58219)  %dorng% {
  
  eps=rnorm(816) # Resample the error term with a normal distribution
  
  # Resample the output
  y_boot=x%*%semipara_full$beta+semipara_full$fitted.values+semipara_full$residuals*eps
  
  # Estimate the model
  optim_full_boot=optim(optim_full2$par,SP_full,sign=-1,y=y_boot, x=x,k=10,method = "BFGS")
  semipara_full_boot = SP_full(optim_full_boot$par,y_boot,x,k=10,opt=FALSE) # Extract best fit model
  lambda_full_boot=predict(semipara_full_boot,type="terms")[,2] # Predicted values of lambda
 
  # Save the results
  results=list(c(semipara_full_boot$alpha,semipara_full_boot$beta),lambda_full_boot,mean(-x%*%semipara_full_boot$beta))
  
}
proc.time() - ptm

stopCluster(cl)


#############################################################################################
####### Wild bootstrap without constraints on alpha and heteroscedasticity correction #######
#############################################################################################

workers <- 10 # Number of splited cores
cl <- parallel::makeCluster(workers)

# register the cluster for using foreach
doParallel::registerDoParallel(cl)

ptm = proc.time()
rep=2 # number of replicates
Bootstrap_full_hete=foreach (i=1:rep,.export=c("gam"),.combine=cbind,.options.RNG = 58219)  %dorng% {
  
  eps=rnorm(816) # Resample the error term with a normal distribution
  
  # Resample the output
  y_boot=x%*%semipara_full_hete$beta+semipara_full_hete$fitted.values+semipara_full_hete$residuals*eps
  
  # Estimate the model
  optim_full_boot=optim(optim_full2$par,SP_full,sign=-1,y=y_boot, x=x,k=10,method = "BFGS")
  semipara_full_boot = SP_full(optim_full_boot$par,y_boot,x,k=10,opt=FALSE) # Extract best fit model
  
  # Estimate the heteroscedasticity
  v2<-abs(semipara_full_boot$residuals)
  hetero_full_boot <-optim(semipara_full_boot$alpha[1:7], hetero.est, x=x,e=v2) 
  hetero.fit_full$par <- hetero.est(hetero_full_boot$par,x,e=v2,opt=FALSE)# Extract best fit model
  w_boot<-1/exp(hetero.fit_full$fitted.values) # Weights that correct for heteroscedasticity
  
  # Estimate the model
  optim_full_boot_hete=optim(optim_full2_hete$par,SP_full,sign=-1,y=y_boot, x=x,w=w_boot,k=10,method = "BFGS")
  semipara_full_boot_hete = SP_full(optim_full_boot_hete$par,y=y_boot,x=x,w=w_boot,k=10,opt=FALSE) # Extract best fit model
  lambda_full_boot_hete=predict(semipara_full_boot_hete,type="terms")[,2]/w_boot # Predicted values of lambda
 
  # Save the results
  results=list(c(semipara_full_boot_hete$alpha,semipara_full_boot_hete$beta),lambda_full_boot_hete,mean(-x%*%semipara_full_boot_hete$beta))
  
}

stopCluster(cl)
proc.time() - ptm

#####################################################################
######### Display the results #######################################
#####################################################################

# Bootstrap replicates of coefficients under unconstrained alpha
alpha_SP=Reduce(cbind,Bootstrap[1,])[1:23,]
beta_SP=Reduce(cbind,Bootstrap[1,])[24:46,]
eta_SP=Reduce(cbind,Bootstrap[3,])

# Bootstrap replicates of coefficients under unconstrained alpha and hetero correction
alpha_SP_hete=Reduce(cbind,Bootstrap_hete[1,])[1:23,]
beta_SP_hete=Reduce(cbind,Bootstrap_hete[1,])[24:46,]
eta_SP_hete=Reduce(cbind,Bootstrap_hete[3,])

# Bootstrap replicates of coefficients under constrained alpha
alpha_SP_full=Reduce(cbind,Bootstrap_full[1,])[1:23,]
beta_SP_full=Reduce(cbind,Bootstrap_full[1,])[24:46,]
eta_SP_full=Reduce(cbind,Bootstrap_full[3,])

# Bootstrap replicates of coefficients under constrained alpha and hetero correction
alpha_SP_full_hete=Reduce(cbind,Bootstrap_full_hete[1,])[1:23,]
beta_SP_full_hete=Reduce(cbind,Bootstrap_full_hete[1,])[24:46,]
eta_SP_full_hete=Reduce(cbind,Bootstrap_full_hete[3,])

# Replicated lambda for the different results
boot_lambda=Reduce(cbind,Bootstrap[2,])
boot_lambda_hete=Reduce(cbind,Bootstrap_hete[2,])
boot_lambda_full=Reduce(cbind,Bootstrap_full[2,])
boot_lambda_full_hete=Reduce(cbind,Bootstrap_full_hete[2,])

# Export results
write.csv(alpha_SP,"alpha_SP.csv")
write.csv(alpha_SP_hete,"alpha_SP_hete.csv")
write.csv(alpha_SP_full,"alpha_SP_full.csv")
write.csv(alpha_SP_full_hete,"alpha_SP_full_hete.csv")
write.csv(beta_SP,"beta_SP.csv")
write.csv(beta_SP_hete,"beta_SP_hete.csv")
write.csv(beta_SP_full,"beta_SP_full.csv")
write.csv(beta_SP_full_hete,"beta_SP_full_hete.csv")
write.csv(eta_SP,"eta_SP.csv")
write.csv(eta_SP_hete,"eta_SP_hete.csv")
write.csv(eta_SP_full,"eta_SP_full.csv")
write.csv(eta_SP_full_hete,"eta_SP_full_hete.csv")
write.csv(boot_lambda,"boot_lambda.csv")
write.csv(boot_lambda_hete,"boot_lambda_hete.csv")
write.csv(boot_lambda_full,"boot_lambda_full.csv")
write.csv(boot_lambda_full_hete,"boot_lambda_full_hete.csv")

# Import results if previously saved
alpha_SP <- read.csv("./Bresults/alpha_SP.csv")
alpha_SP_hete <- read.csv("./Bresults/alpha_SP_hete.csv")
alpha_SP_full <- read.csv("./Bresults/alpha_SP_full.csv")
alpha_SP_full_hete <- read.csv("./Bresults/alpha_SP_full_hete.csv")
alpha_SP$X <- NULL
alpha_SP_hete$X <- NULL
alpha_SP_full$X <- NULL
alpha_SP_full_hete$X <- NULL
beta_SP <- read.csv("./Bresults/beta_SP.csv")
beta_SP_hete <- read.csv("./Bresults/beta_SP_hete.csv")
beta_SP_full <- read.csv("./Bresults/beta_SP_full.csv")
beta_SP_full_hete <- read.csv("./Bresults/beta_SP_full_hete.csv")
beta_SP$X <- NULL
beta_SP_hete$X <- NULL
beta_SP_full$X <- NULL
beta_SP_full_hete$X <- NULL
eta_SP <- read.csv("./Bresults/eta_SP.csv")
eta_SP_hete <- read.csv("./Bresults/eta_SP_hete.csv")
eta_SP_full <- read.csv("./Bresults/eta_SP_full.csv")
eta_SP_full_hete <- read.csv("./Bresults/eta_SP_full_hete.csv")
eta_SP$X <- NULL
eta_SP_hete$X <- NULL
eta_SP_full$X <- NULL
eta_SP_full_hete$X <- NULL
boot_lambda <- read.csv("./Bresults/boot_lambda.csv")
boot_lambda_hete <- read.csv("./Bresults/boot_lambda_hete.csv")
boot_lambda_full <- read.csv("./Bresults/boot_lambda_full.csv")
boot_lambda_full_hete <- read.csv("./Bresults/boot_lambda_full_hete.csv")

#####################################
# Bootstrap test for orthogonality ##
#####################################
N<-nrow(data)
# Compute the replicated statistic alpha'beta
boot_test_SP=sqrt(N)*diag(t(alpha_SP)%*%as.matrix(beta_SP))
boot_test_SP_hete=sqrt(N)*diag(t(alpha_SP_hete)%*%as.matrix(beta_SP_hete))
boot_test_SP_full=sqrt(N)*diag(t(alpha_SP_full)%*%as.matrix(beta_SP_full))
boot_test_SP_full_hete=sqrt(N)*diag(t(alpha_SP_full_hete)%*%as.matrix(beta_SP_full_hete))

# Standars error of the statistic
v_test_SP=var(boot_test_SP)
v_test_SP_hete=var(boot_test_SP_hete)
v_test_SP_full=var(boot_test_SP_full)
v_test_SP_full_hete=var(boot_test_SP_full_hete)

# Compute the statistic alpha'beta
test_SP=N*(semipara$alpha%*%semipara$beta)^2/v_test_SP
test_SP_hete=N*(semipara_hete$alpha%*%semipara_hete$beta)^2/v_test_SP_hete
test_SP_full=N*(semipara_full$alpha%*%semipara_full$beta)^2/v_test_SP_full
test_SP_full_hete=N*(semipara_full_hete$alpha%*%semipara_full_hete$beta)^2/v_test_SP_full_hete


#####################
# Tables of results #
#####################

# Reported Coefficient Estimates and standard error from Orea at al. (2015)
Estimate_exp=c(5.042,0.236,-0.114,0.799,-1.469,0.348,0.080,0.005,-6.148,1.446,6.148,-12.187,0,0,0)
SE_exp=c(0.018,0.046,0.030,0.047,0.088,0.013,0.008,0.001,2.034,0.769,2.034,3.165,0,0,0)
Estimate_logis=c(5.043,0.238,-0.117,0.797,-1.480,0.347,0.080,.005,-7.014,1.577,7.014,-14.283,0,0,0)
SE_logis=c(0.018,0.046,0.030,0.047,0.086,0.013,0.008,0.001,2.242,0.862,2.242,3.640,0,0,0)

# Coefficient Estimates and standard error under unconstrained alpha
Estimate_SP=c(mean(-x%*%semipara$beta),semipara$beta[1:7],semipara$alpha[1:7])
SE_SP=c(sqrt(apply(eta_SP,1,var)),sqrt(apply(beta_SP,1,var))[1:7],sqrt(apply(alpha_SP,1,var))[1:7])

# Coefficient Estimates and standard error under unconstrained alpha and hetero correction
Estimate_SP_hete=c(mean(-x%*%semipara_hete$beta),semipara_hete$beta[1:7],semipara_hete$alpha[1:7])
SE_SP_hete=c(sqrt(apply(eta_SP_hete,1,var)),sqrt(apply(beta_SP_hete,1,var))[1:7],sqrt(apply(alpha_SP_hete,1,var))[1:7])

# Coefficient Estimates and standard error  under constrained alpha
Estimate_SP_full=c(mean(-x%*%semipara_full$beta),semipara_full$beta[1:7],semipara_full$alpha[1:7])
SE_SP_full=c(sqrt(apply(eta_SP_full,1,var)),sqrt(apply(beta_SP_full,1,var))[1:7],sqrt(apply(alpha_SP_full,1,var))[1:7])

# Coefficient Estimates and standard error  under constrained alpha and hetero correction
Estimate_SP_full_hete=c(mean(-x%*%semipara_full_hete$beta),semipara_full_hete$beta[1:7],semipara_full_hete$alpha[1:7])
SE_SP_full_hete=c(sqrt(apply(eta_SP_full_hete,1,var)),sqrt(apply(beta_SP_full_hete,1,var))[1:7],sqrt(apply(alpha_SP_full_hete,1,var))[1:7])

###### Significance of the estimated coefficients of the constrained model

# Significance at 1%
apply(beta_SP,1,quantile,c(0.005,0.995))[,1:7] #All significant
apply(alpha_SP,1,quantile,c(0.005,0.995))[,1:7] #All significant
apply(eta_SP,1,quantile,c(0.005,0.995)) #Not significant

# Significance at 5%
#apply(beta_SP,1,quantile,c(0.025,0.975))[,1:7] #All significant
#apply(alpha_SP,1,quantile,c(0.025,0.975))[,1:7] #All significant
apply(eta_SP,1,quantile,c(0.025,0.975))  #Significant

###### Significance of the estimated coefficients of the constrained model with hete correction

# Significance at 1%
apply(beta_SP_hete,1,quantile,c(0.005,0.995))[,1:7] #All significant
apply(alpha_SP_hete,1,quantile,c(0.005,0.995))[,1:7] #All significant
apply(eta_SP_hete,1,quantile,c(0.005,0.995)) #Not significant

# Significance at 5%
#apply(beta_SP_hete,1,quantile,c(0.025,0.975))[,1:7] #All significant
#apply(alpha_SP_hete,1,quantile,c(0.025,0.975))[,1:7] #All significant
apply(eta_SP_hete,1,quantile,c(0.025,0.975)) #Not significant

# Significance at 10%
#apply(beta_SP_hete,1,quantile,c(0.05,0.95))[,1:7] #All significant
#apply(alpha_SP_hete,1,quantile,c(0.05,0.95))[,1:7] #All significant
apply(eta_SP_hete,1,quantile,c(0.05,0.95)) #Not significant

###### Significance of the estimated coefficients of the unconstrained model

# Significance at 1%
apply(beta_SP_full,1,quantile,c(0.005,0.995))[,1:7] #All except ln(HDD)
apply(alpha_SP_full,1,quantile,c(0.005,0.995))[,1:7] #All except ln(P)
apply(eta_SP_full,1,quantile,c(0.005,0.995)) #Significant

# Significance at 5%
apply(beta_SP_full,1,quantile,c(0.025,0.975))[,1:7] #All significant
apply(alpha_SP_full,1,quantile,c(0.025,0.975))[,1:7] #All except ln(P)
apply(eta_SP_full,1,quantile,c(0.025,0.975)) #Significant

# Significance at 10%
apply(beta_SP_full,1,quantile,c(0.05,0.95))[,1:7] #All significant
apply(alpha_SP_full,1,quantile,c(0.05,0.95))[,1:7] #All except ln(P)
apply(eta_SP_full,1,quantile,c(0.05,0.95)) #Significant

###### Significance of the estimated coefficients of the unconstrained model with hete correction
# Significance at 1%
apply(beta_SP_full_hete,1,quantile,c(0.005,0.995))[,1:7] #All except ln(HDD)
apply(alpha_SP_full_hete,1,quantile,c(0.005,0.995))[,1:7] #All except ln(P)
apply(eta_SP_full_hete,1,quantile,c(0.005,0.995)) #Not significant

# Significance at 5%
apply(beta_SP_full_hete,1,quantile,c(0.025,0.975))[,1:7] #All significant
apply(alpha_SP_full_hete,1,quantile,c(0.025,0.975))[,1:7] #All except ln(P)
apply(eta_SP_full_hete,1,quantile,c(0.025,0.975)) #Not significant

# Significance at 10%
apply(beta_SP_full_hete,1,quantile,c(0.05,0.95))[,1:7] #All significant
apply(alpha_SP_full_hete,1,quantile,c(0.05,0.95))[,1:7] #All except ln(P)
apply(eta_SP_full_hete,1,quantile,c(0.05,0.95)) #Not significant

# build table for coefficients
Result_coefficients = data.table("nm_var" = c("Intercept","ln(Y)","ln(P)","ln(POP)","ln(AHS)","ln(HDD)","ln(CDD)","SDH","ln(Y)","ln(P)",
                                              "ln(POP)","ln(AHS)","ln(HDD)","ln(CDD)","SDH"),
                                 "Est."=c(Estimate_logis),"SE"=c(SE_logis),"Est."=c(Estimate_exp),"SE"=c(SE_exp),"Est."=Estimate_SP,"SE"=SE_SP 
                                 ,"Est."=Estimate_SP_hete,"SE"=SE_SP_hete,"Est."=Estimate_SP_full,"SE"=SE_SP_full,
                                 "Est."=Estimate_SP_full_hete,"SE"=SE_SP_full_hete)

# build table for orthogonality test
Result_tests = data.table("method" = c("PS(r)","WPS(r)","PS(u)","WPS(u)")
                          ,"Est."=c(test_SP,test_SP_hete,test_SP_full,test_SP_full_hete))

#Display Tables
print(Result_coefficients,digits = 3)
print(t(Result_tests),digits = 3)

# Tex format
xtable(Result_coefficients, type = "latex",digits = 3)
xtable(t(Result_tests), type = "latex",digits = 3)

##################
##### Graphs #####
##################

########## Parametric with exponential link #######

y_SC=(y-Result_ML_exp$beta0- x%*%Result_ML_exp$beta) # Estimated lambda + residuals
SC2=sqrt(2/pi)*(1-(exp(-Result_ML_exp$alpha0+x%*%Result_ML_exp$alpha)-1)/exp(-Result_ML_exp$alpha0+x%*%Result_ML_exp$alpha))
a=-Result_ML_exp$alpha0+x%*%Result_ML_exp$alpha
graph_SC=data.frame(a,SC2,y_SC)
ggplot(data=graph_SC,aes(x=a,y=SC2))+
  geom_line(size=2)+theme_bw(base_family="Helvetica")+theme_classic(base_size=15)+
  geom_point(aes(x=a,y=y_SC),size = 1)+
  theme(axis.text.y = element_text(angle = 90, hjust = 0.8,size=30), 
        axis.text.x = element_text(size=30),
        axis.title.x = element_text(vjust = - 0.5))+
  xlab(TeX("$X' \\widehat{\\alpha}_{exponential}$"))+
  ylab(TeX("$Y-X' \\widehat{\\beta}_{exponential}$"))+
  theme(axis.title=element_text(size=40,face="bold"))


########## Parametric with logistic link #######

y_PA=(y-Result_ML_logis$beta0- x%*%Result_ML_logis$beta) # Estimated lambda + residuals
PA2=sqrt(2/pi)*(1-exp(-Result_ML_logis$alpha0+x%*%Result_ML_logis$alpha)/(1+exp(-Result_ML_logis$alpha0+x%*%Result_ML_logis$alpha)))
a=-Result_ML_logis$alpha0+x%*%Result_ML_logis$alpha
graph_PA=data.frame(a,PA2,y_PA)
ggplot(data=graph_PA,aes(x=a,y=PA2))+
  geom_line(size=2)+theme_bw(base_family="Helvetica")+theme_classic(base_size=15)+
  geom_point(aes(x=a,y=y_PA),size = 1)+
  theme(axis.text.y = element_text(angle = 90, hjust = 0.8,size=30), 
        axis.text.x = element_text(size=30),
        axis.title.x = element_text(vjust = - 0.5))+
  xlab(TeX("$X' \\widehat{\\alpha}_{logistic}$"))+
  ylab(TeX("$Y-X' \\widehat{\\beta}_{logistic}$"))+
  theme(axis.title=element_text(size=40,face="bold"))



########### Semi-parametric with constrained alpha

PS2=lambda-min(lambda)  # Estimated lambda with minimum at 0
x_PS=x%*%semipara$alpha
y_PS=(y- x%*%semipara$beta)-min(lambda) # Estimated lambda + residuals
upr=apply(boot_lambda,1,quantile,0.975)-min(lambda) # Upper bound for 95% CI
lwr=apply(boot_lambda,1,quantile,0.025)-min(lambda) # Lower bound for 95% CI

# Build graph
graph_PS=data.frame(x_PS,PS2,y_PS,lwr,upr)
ggplot(data=graph_PS,aes(x=x_PS,y=PS2))+
  geom_ribbon(aes(x=x_PS,ymin=lwr, ymax=upr), show.legend=FALSE, fill="gray70", alpha=0.5)+
  geom_line(size=2)+theme_bw(base_family="Helvetica")+theme_classic(base_size=15)+
  geom_point(aes(x=x_PS,y=y_PS),size = 1)+
  theme(axis.text.y = element_text(angle = 90, hjust = 0.8,size=30), 
        axis.text.x = element_text(size=30),
        axis.title.x = element_text(vjust = - 0.5))+
  xlab(TeX("$X' \\widehat{\\alpha}_{PS}$"))+
  ylab(TeX("$Y-X' \\widehat{\\beta}_{PS}$"))+
  theme(axis.title=element_text(size=40,face="bold"))+
  xlim(-0.20, 0.227)+ylim(-0.134, 0.44)

# Here are computed to bounds of the axes
#min(graph_PS$x_PS)
#min(graph_PS_hete$x_PS_hete)
#max(graph_PS$x_PS)
#max(graph_PS_hete$x_PS_hete)
#min(graph_PS$y_PS)
#min(graph_PS_hete$y_PS_hete)
#max(graph_PS$y_PS)
#max(graph_PS_hete$y_PS_hete)

########### Semi-parametric with constrained alpha with heteroscedasticity correction

PS2_hete=lambda_hete-min(lambda_hete)  # Estimated lambda with minimum at 0
x_PS_hete=x%*%semipara_hete$alpha
y_PS_hete=(y- x%*%semipara_hete$beta)-min(lambda_hete) # Estimated lambda + residuals
upr_hete=apply(boot_lambda_hete,1,quantile,0.975)-min(lambda_hete) # Upper bound for 95% CI
lwr_hete=apply(boot_lambda_hete,1,quantile,0.025)-min(lambda_hete) # Lower bound for 95% CI

# Build graph
graph_PS_hete=data.frame(x_PS_hete,PS2_hete,y_PS_hete,lwr_hete,upr_hete)
ggplot(data=graph_PS_hete,aes(x=x_PS_hete,y=PS2_hete))+
  geom_ribbon(aes(x=x_PS_hete,ymin=lwr_hete, ymax=upr_hete), show.legend=FALSE, fill="gray70", alpha=0.5)+
  geom_line(size=2)+theme_bw(base_family="Helvetica")+theme_classic(base_size=15)+
  geom_point(aes(x=x_PS_hete,y=y_PS_hete),size = 1)+
  theme(axis.text.y = element_text(angle = 90, hjust = 0.8,size=30), 
        axis.text.x = element_text(size=30),
        axis.title.x = element_text(vjust = - 0.5))+
  xlab(TeX("$X' \\widehat{\\alpha}_{WPS}$"))+
  ylab(TeX("$Y-X' \\widehat{\\beta}_{WPS}$"))+
  theme(axis.title=element_text(size=40,face="bold"))+
  xlim(-0.20, 0.227)+ylim(-0.134, 0.44)


########### Semi-parametric with unconstrained alpha

PS2_full=lambda_full-min(lambda_full) # Estimated lambda with minimum at 0
x_PS_full=x%*%semipara_full$alpha
y_PS_full=(y- x%*%semipara_full$beta)-min(lambda_full) # Estimated lambda + residuals
upr_full=apply(boot_lambda_full,1,quantile,0.975)-min(lambda_full) # Upper bound for 95% CI
lwr_full=apply(boot_lambda_full,1,quantile,0.025)-min(lambda_full) # Lower bound for 95% CI

# Build graph
graph_PS_full=data.frame(x_PS_full,PS2_full,y_PS_full,lwr_full,upr_full)
ggplot(data=graph_PS_full,aes(x=x_PS_full,y=PS2_full))+
  geom_ribbon(aes(x=x_PS_full,ymin=lwr_full, ymax=upr_full), show.legend=FALSE, fill="gray70", alpha=0.5)+
  geom_line(size=2)+theme_bw(base_family="Helvetica")+theme_classic(base_size=15)+
  geom_point(aes(x=x_PS_full,y=y_PS_full),size = 1)+
  theme(axis.text.y = element_text(angle = 90, hjust = 0.8,size=30), 
        axis.text.x = element_text(size=30),
        axis.title.x = element_text(vjust = - 0.5))+
  xlab(TeX("$X' \\widehat{\\alpha}$"))+
  ylab(TeX("$Y-X' \\widehat{\\beta}$"))+
  theme(axis.title=element_text(size=40,face="bold"))+
  xlim(-1.391,0.771)+ylim(-0.153, 1.263)

#Here are computed to bounds of the axes
#min(graph_PS_full$x_PS_full)
#min(graph_PS_full_hete$x_PS_full_hete)
#max(graph_PS_full$x_PS_full)
#max(graph_PS_full_hete$x_PS_full_hete)
#min(graph_PS_full$lwr_full)
#min(graph_PS_full_hete$lwr_full_hete)
#max(graph_PS_full$upr_full)
#max(graph_PS_full_hete$upr_full_hete)

########### Semi-parametric with unconstrained alpha with heteroscedasticity correction

PS2_full_hete=lambda_full_hete-min(lambda_full_hete) # Estimated lambda with minimum at 0
x_PS_full_hete=x%*%semipara_full_hete$alpha
y_PS_full_hete=(y- x%*%semipara_full_hete$beta)-min(lambda_full_hete) # Estimated lambda + residuals
upr_full_hete=apply(boot_lambda_full_hete,1,quantile,0.975)-min(lambda_full_hete) # Upper bound for 95% CI
lwr_full_hete=apply(boot_lambda_full_hete,1,quantile,0.025)-min(lambda_full_hete) # Lower bound for 95% CI

# Build graph
graph_PS_full_hete=data.frame(x_PS_full_hete,PS2_full_hete,y_PS_full_hete,lwr_full_hete,upr_full_hete)
ggplot(data=graph_PS_full_hete,aes(x=x_PS_full_hete,y=PS2_full_hete))+
  geom_ribbon(aes(x=x_PS_full_hete,ymin=lwr_full_hete, ymax=upr_full_hete), show.legend=FALSE, fill="gray70", alpha=0.5)+
  geom_line(size=2)+theme_bw(base_family="Helvetica")+theme_classic(base_size=15)+
  geom_point(aes(x=x_PS_full_hete,y=y_PS_full_hete),size = 1)+
  theme(axis.text.y = element_text(angle = 90, hjust = 0.8,size=30), 
        axis.text.x = element_text(size=30),
        axis.title.x = element_text(vjust = - 0.5))+
  xlab(TeX("$X' \\widehat{\\alpha}_{W}$"))+
  ylab(TeX("$Y-X' \\widehat{\\beta}_{W}$"))+
  theme(axis.title=element_text(size=40,face="bold"))+
  xlim(-1.391,0.771)+ylim(-0.153, 1.263)




# plot of rebounds
ggplot(data=graph_SC,aes(x=a,y=1-SC2))+geom_line(size=2,linetype='dotted')+
  geom_line(data=graph_PA,aes(x=a,y=1-PA2),size=2)+
  theme_classic(base_size=15)+
  geom_hline(yintercept = 1.1)+
  geom_hline(yintercept = 1)+
  geom_hline(yintercept = .9)+
  geom_hline(yintercept = .8)+
  geom_hline(yintercept = .7)+
  geom_hline(yintercept = .6)+
  geom_hline(yintercept = .5)+
  geom_hline(yintercept = .4)+
  geom_hline(yintercept = .3)+
  geom_hline(yintercept = .2)+
  geom_hline(yintercept = .1)+
  ylab("Rebound")+xlab(TeX("$X'\\hat{\\alpha}$"))+ theme(axis.title=element_text(size=40,face="bold"))+
  ylim(0, 1.1)+
  scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),breaks = seq(0,1.1,by=.1))

ggplot(data=graph_SC,aes(x=a,y=1-SC2))+
  geom_line(data=graph_PS,aes(x=x_PS,y=1-PS2),size=2)+
  theme_classic(base_size=15)+
  geom_hline(yintercept = 1.1)+
  geom_hline(yintercept = 1)+
  geom_hline(yintercept = .9)+
  geom_hline(yintercept = .8)+
  geom_hline(yintercept = .7)+
  geom_hline(yintercept = .6)+
  geom_hline(yintercept = .5)+
  geom_hline(yintercept = .4)+
  geom_hline(yintercept = .3)+
  geom_hline(yintercept = .2)+
  geom_hline(yintercept = .1)+
  ylab("Rebound")+ xlab(TeX("$X'\\hat{\\alpha}$"))+ theme(axis.title=element_text(size=40,face="bold"))+
  ylim(0, 1.1)+
  scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),breaks = seq(0,1.1,by=.1))


ggplot(data=graph_SC,aes(x=a,y=1-SC2))+
  geom_line(data=graph_PS_full,aes(x=x_PS_full,y=1-PS2_full),size=2)+
  theme_classic(base_size=15)+
  geom_hline(yintercept = 1.1)+
  geom_hline(yintercept = 1)+
  geom_hline(yintercept = .9)+
  geom_hline(yintercept = .8)+
  geom_hline(yintercept = .7)+
  geom_hline(yintercept = .6)+
  geom_hline(yintercept = .5)+
  geom_hline(yintercept = .4)+
  geom_hline(yintercept = .3)+
  geom_hline(yintercept = .2)+
  geom_hline(yintercept = .1)+
  ylab("Rebound")+ xlab(TeX("$X'\\hat{\\alpha}$"))+ theme(axis.title=element_text(size=40,face="bold"))+
  ylim(0, 1.1)+
  scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),breaks = seq(0,1.1,by=.1))

####################
# comparisons
###############

x_SC<--Result_ML_exp$alpha0+x%*%Result_ML_exp$alpha
x_PA<--Result_ML_logis$alpha0+x%*%Result_ML_logis$alpha

temp<-data.frame(y_PA,y_SC,y_PS,y_PS_full)
cor(temp)

temp2<-data.frame(x_PA,x_SC,x_PS,x_PS_full)
cor(temp2)

temp<-data.frame(y_PA,y_SC,y_PS,y_PS_full,x_PA,x_SC,x_PS,x_PS_full)



############################################################
# Figure 2
############################################################
library(ggplot2)
library(scam)


t<-scam(y_PS~s(x_PS,bs="mdcx",k=10))
yhat_PS<-predict(t,type="response")

t<-gam(y_PS~s(x_PS,bs="cr",k=10),method="REML")
yhat_PS_u<-predict(t,type="response")

qplot(x_PS,y_PS)+
  geom_line(aes(x=x_PS,y=yhat_PS),size=2,linetype="dotdash")+
  geom_line(aes(x=x_PS,y=yhat_PS_u),size=2)

temp2<-data.frame(x_PS, y_PS, yhat_PS, yhat_PS_u)

ggplot(temp2)+aes(x=x_PS,y=y_PS)+geom_point(size=.1)+
  geom_line(aes(x=x_PS,y=yhat_PS),size=1,linetype="dashed")+
  geom_line(aes(x=x_PS,y=yhat_PS_u),size=1)+
  theme_classic(base_size=15)+
  ylab(TeX("$Y-X' \\widehat{\\beta}_{PS}$"))+
  xlab(TeX("$X'\\hat{\\alpha}_{PS}$"))+ 
  theme(axis.title=element_text(size=40,face="bold"))
  


t<-scam(y_PS_full~s(x_PS_full,bs="mpd",k=10))
yhat_PS_full<-predict(t,type="response")

t<-gam(y_PS_full~s(x_PS_full,bs="cr",k=10),method="REML")
yhat_PS_full_u<-predict(t,type="response")
                        
t<-gam(y_PS_full~s(x_PS_full,bs="cr",k=10),sp=100)
yhat_PS_full_uu<-predict(t,type="response")
temp2<-data.frame(x_PS_full, y_PS_full, yhat_PS_full, yhat_PS_full_u, yhat_PS_full_uu)


ggplot(temp2)+aes(x=x_PS_full,y=y_PS_full)+geom_point(size=.1)+
  geom_line(aes(x=x_PS_full,y=yhat_PS_full),size=1,linetype="dashed")+
  geom_line(aes(x=x_PS_full,y=yhat_PS_full_u),size=1)+
  geom_line(aes(x=x_PS_full,y=yhat_PS_full_uu),size=1,linetype="dotdash")+
  theme_classic(base_size=15)+
  ylab(TeX("$Y-X' \\widehat{\\beta}$"))+
  xlab(TeX("$X'\\hat{\\alpha}$"))+ 
  theme(axis.title=element_text(size=40,face="bold"))

