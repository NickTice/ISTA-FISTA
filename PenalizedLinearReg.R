# Penalized linear regression with LASSO, SCAD, and MCP penalties

penalizedLin=function(
  X,
  y,
  penalty,
  lambda,
  gamma=NA,
  tol=1e-8,
  max.iter=1000
){
  gamma=ifelse(is.na(gamma) & penalty=="SCAD",3.7,gamma )
  gamma=ifelse(is.na(gamma) & penalty=="MCP",3,gamma )
  # Number of observations
  n = length(y)
  
  # Number of features
  d = dim(X)[2]
  
  ## Center y, but do **NOT** scale y
  y.mean = mean(y) ## Need to keep track of the mean of y for later
  y = y-y.mean # Subtract y.mean from every entry of y
  
  ## Center X.
  X.means = colMeans(X) ## Need to keep track of the column means for later
  X = sweep(X,2,X.means, FUN='-') # Subtract every entry in each column of X by the corresponding column mean
  
  ## Take note of the column population standard deviations, so we can transform back later
  #
  ## Note that the sd function in R gives the **sample** standard deviation, but we really want
  ## the population standard deviation, so we multiply by sqrt((n-1)/n
  sigma.X = sqrt((n-1)/n) * apply(X,2,sd)
  
  ## Divide each jth column of X by the jth element in sigma.X
  X = sweep(X, 2, sigma.X, FUN = '/') 

# soft thresholding operator
 S<-function(zj,lam){
  if (zj>lam) return(zj-lam)
  if (abs(zj)<=lam) return(0)
  if (zj<(-1)*lam) return(zj+lam)
}

  fscad<-function(zj, lambda){
    if (abs(zj)<=2*lambda) return(S(zj,lambda))
    if (2*lambda<abs(zj) && abs(zj)<=gamma*lambda) return((gamma-1)/(gamma-2)*S(zj,(gamma*lambda)/(gamma-1)))
    if (abs(zj)>gamma*lambda) return(zj)
  }
  
  fmcp<-function(zj, lambda){
    if (abs(zj)<=gamma*lambda) return(gamma/(gamma-1)*S(zj,lambda))
    if (abs(zj)>gamma*lambda) return(zj)
  }
  
  betas=c()
  beta0s=c()
  

  
  for (lam in lambda) {
    
    if (penalty=="lasso"){
      
      z=rep(0,d)
      beta1=rep(0,d)
      beta2=rep(1,d)
      
      
      r<- y-X%*%beta1
      i<-1
      while (i<=max.iter && max(abs(beta1-beta2))>tol) {
        beta1<-beta2
        for (j in 1:d) {
          z[j]<-beta1[j]+(1/n)*t(X[,j])%*%r
          beta2[j]<-S(z[j], lam)
          r=r-(beta2[j]-beta1[j])*X[,j]
        }
        i<-i+1
      }
      
    } 
    
    if (penalty=="SCAD"){
      
      z=rep(0,d)
      beta1=rep(0,d)
      beta2=rep(1,d)
      
      ## SCAD 
      
      r<- y-X%*%beta1
      i<-1
      while (i<=max.iter && max(abs(beta1-beta2))>tol) {
        beta1<-beta2
        for (j in 1:d) {
          z[j]<-beta1[j]+(1/n)*t(X[,j])%*%r
          beta2[j]<-fscad(z[j], lam)
          r=r-(beta2[j]-beta1[j])*X[,j]
        }
        i<-i+1
      }
      
    } 
    
    
    ## MCP
    
    if (penalty=="MCP"){
      
      z=rep(0,d)
      beta1=rep(0,d)
      beta2=rep(1,d)
      
      ## MCP
      
      r<- y-X%*%beta1
      i<-1
      while (i<=max.iter && max(abs(beta1-beta2))>tol) {
        beta1<-beta2
        for (j in 1:d) {
          z[j]<-beta1[j]+(1/n)*t(X[,j])%*%r
          beta2[j]<-fmcp(z[j], lam)
          r=r-(beta2[j]-beta1[j])*X[,j]
        }
        i<-i+1
      }
      
    } 
    
    ## Now we need to transform beta2 BACK to the original scale and get the intercept
    
    ## beta on the original scale of the features
    beta = beta2 / sigma.X
    beta0 = y.mean - sum(beta*X.means)
    
    betas=cbind(betas, beta)
    beta0s=cbind(beta0s, beta0)
    
  }
  
  colnames(betas)=lambda
  colnames(beta0s)=lambda
  
  
  ## Return beta and beta0
  # Return list   
  return(list(beta = betas,
              beta0 = beta0s))
  
}
