# accel=FALSE for ISTA and accel=TRUE for FISTA

ISTA=function(X,
                   y,
                   t,
                   lambda,
                   accel,
                   tol=1e-5,
                   max.iter=5000){
  
  # Soft thresholding
  S<-function(u,lam){
    if (u>lam) return(u-lam)
    if ((u <= lam) & (u>= (-1)*lam)) return(0)
    if (u<(-1)*lam) return(u+lam)
  }
  
  n = length(y)
  
  d = dim(X)[2]
  
  # standarization
  y.mean = mean(y)
  y = y-y.mean
  
  X.means = colMeans(X) 
  X = sweep(X,2,X.means, FUN='-')
  
  sigma.X = sqrt((n-1)/n) * apply(X,2,sd)
  
  X = sweep(X, 2, sigma.X, FUN = '/') 
  
  if (accel==FALSE){
    # ISTA
    beta_k=rep(1,d)
    beta_k1=rep(0,d)
    i=0
    fs.hat=c()
    while (i<=max.iter && max(abs(beta_k-beta_k1))>tol) {
      
      beta_k=beta_k1
      
      # applying the algorithm
      beta_k1=sapply(beta_k-t*1/n*t(X) %*% (X %*% beta_k - y), S, lam=lambda*t)
      
      # calculate fs hats
      fs=1/(2*n)*(norm(y-X%*%beta_k1, "2"))^2+lambda*norm(as.matrix(beta_k1), "O")
      
      fs.hat=c(fs.hat, fs)
      
      i=i+1
    }
    betas=beta_k1
  }
  
  if (accel==TRUE){
    # FISTA
    beta_x0=rep(1,d)
    beta_x1=rep(1,d)
    beta_x2=rep(0,d)
    i=1
    k=2
    fs.hat=c()
    
    
    while (i<=max.iter && max(abs(beta_x2-beta_x1))>tol) {
      
      beta_x0=beta_x1
      beta_x1=beta_x2
      
      # momentum
      v=beta_x1+(k-2)/(k+1)*(beta_x1-beta_x0)
      
      # apply algorithm
      beta_x2=sapply((v-t*1/n*t(X) %*% (X %*% v - y)), S, lam=lambda*t)
      
      k=k+1
      
      fs=1/(2*n)*(norm(y-X%*%beta_x2, "2"))^2+lambda*norm(as.matrix(beta_x2), "O")
      
      fs.hat=c(fs.hat, fs)
      
      
      i=i+1
      
    }
    
    betas=beta_x2
  }
  
  # rescale
  beta = betas / sigma.X
  
  # calculate intercept
  beta0 = y.mean - sum(beta*X.means)
  
  return(list(beta=beta, beta0=beta0, fs.hat=fs.hat))
}
