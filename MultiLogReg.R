# Multinomial Logistic Regression with Gradient Descent Algorithm

MultiLogReg=function(X,
                  y,
                  t,
                  lambda,
                  tol=1e-3,
                  max.iter=1000){
  n = length(y)
  
  d = dim(X)[2]
  
  
  X.means = colMeans(X) 
  X = sweep(X,2,X.means, FUN='-')
  
  sigma.X = sqrt((n-1)/n) * apply(X,2,sd)
  
  # If column standard deviation is 0, replace it with 1 so we don't have Inf's
  sigma.X = replace(sigma.X, sigma.X==0, 1)
  
  X = sweep(X, 2, sigma.X, FUN = '/') 
  
  X=cbind(rep(1,n),X)
  
  K=unique(y)
  K=sort(K)
  
  X[is.na(X)]=0
  
  # Matrix of betas. Each column is a different K
  betas=matrix(rep(0, length(K)*(d+1)), nrow = (d+1), ncol = length(K))
  
  it=1
  convergence=1
  while (it<=max.iter && convergence>tol) {
    
    convergence=c()
    for (i in 1:length(K)) {
      beta_i=betas[,i]
      
      k=K[i]
      
      # different terms of the algorithm from part c
      term1=lambda*beta_i
      
      term2=1/n*colSums(ifelse(y==k,1,0)*X, na.rm = TRUE)
      
      y_j=as.matrix(X) %*% beta_i
      
      sm.sum=rep(0,length(y_j))
      
      for (z in 1:length(y_j)) {
        sm.sum[z]=sum(exp(y_j-y_j[z]))
      }
      
      sm=1/(1+sm.sum)
      
      term3=1/n*colSums(sm %*% as.matrix(X), na.rm = TRUE)
      
      
      beta_new=beta_i-t*(term1 - term2 + as.numeric(term3))
      
      betas[,i]=beta_new
      convergence=c(convergence, max(abs(beta_i-beta_new), na.rm = TRUE))
    }
    
    convergence=max(convergence)
    it=it+1
    
    # print iteration number and max convergence
    print(it)
    print(convergence)
  }
  
  # scale back
  sigma.X=ifelse(sigma.X==0, 1, sigma.X)  
  
  for (i in 1:ncol(betas)) {
    betas[-1,i]=betas[-1,i]/sigma.X
  }
  
  return(betas)
  
}