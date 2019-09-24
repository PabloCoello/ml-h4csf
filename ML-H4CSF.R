if(!require(csn)){install.packages("csn");library(csn)} 
if(!require(stats4)){install.packages("stats4");library(stats4)} 
if(!require(bbmle)){install.packages("bbmle");library(bbmle)} 
if(!require(matrixcalc)){install.packages("matrixcalc");library(matrixcalc)} 


Y=c(2,4,3,5,6)
X=c(1,3,6,2,5)
Z=c(3,5,8,5,6)

log_likelihood = function(beta1, sigmaesq, gamma0, gamma1){
  t=length(Y)
  get_p_matrix = function(t){
    p = -diag(t)
    p = p[-nrow(p),]
    aux = diag(t)
    aux = aux[-1,]
    p = p+aux
    return(p)
  }
  
  P = get_p_matrix(t)
  y = P%*%Y
  x = P%*%X
  R = y - x*beta1
  
  var_array = array(dim=t)
  for (i in 1:t){
    var_array[i]=exp(gamma0+gamma1*Z[i])
  }
  
  
  mu = replicate(t-1, 0)
  sigma = as.matrix(P%*%(sigmaesq*diag(t) + diag(var_array))%*%t(P))
  rownames(sigma) <- colnames(sigma)
  gamma = -diag(var_array)%*%t(P)%*%solve(sigma)
  nu = replicate(t, 0)
  delta = diag(var_array)-diag(var_array)%*%t(P)%*%solve(sigma)%*%P%*%diag(var_array)
  
  R = dcsn(x=array(R), mu, sigma, gamma, nu, delta)
  -sum(log(R))
}

mle2(log_likelihood, start = list(beta1 = 1, sigmaesq = 1, gamma0 = 1, gamma1=1),
     method = "L-BFGS-B",
     trace = TRUE, 
     lower = c(beta1 = -Inf, sigmaesq = 0.001, gamma0 = -Inf, gamma1= -Inf),
     upper = c(beta1 = Inf, sigmaesq = 10, gamma0 = Inf, gamma1= Inf))




