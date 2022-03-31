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
library(lmboot)
library(foreach)
library(doRNG)
library(snow)
library(doParallel)
library(data.table)


###########################################
####### Estimation functions ##############
###########################################

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

load("data_for_super_computer.RData")
#adress="data.txt"
#data = read.csv(adress, sep="")


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


#######################################################
####### Wild bootstrap with constraints on alpha ######
#######################################################

###### Make the core works in parallel ########
workers <- 28
cl <- parallel::makeCluster(workers)
# register the cluster for using foreach
doParallel::registerDoParallel(cl)


rep=1000
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




##########################################################################################
####### Wild bootstrap with constraints on alpha and heteroscedasticity correction #######
##########################################################################################


 # number of replicates
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






##########################################################
####### Wild bootstrap without constraints on alpha ######
##########################################################


# number of replicates
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





#############################################################################################
####### Wild bootstrap without constraints on alpha and heteroscedasticity correction #######
#############################################################################################


 # number of replicates
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



