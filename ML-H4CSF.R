if(!require(csn)){install.packages("csn");library(csn)} 
if(!require(stats4)){install.packages("stats4");library(stats4)} 
if(!require(bbmle)){install.packages("bbmle");library(bbmle)} 
if(!require(matrixcalc)){install.packages("matrixcalc");library(matrixcalc)} 
if(!require(xlsx)){install.packages("xlsx");library(xlsx)}
if(!require(readxl)){install.packages("readxl");library(readxl)} 
if(!require(EnvStats)){install.packages("EnvStats");library(EnvStats)} 
if(!require(plm)){install.packages("plm");library(plm)} 
if(!require(npsf)){install.packages("npsf");library(npsf)} 

setwd("~/GitHub/data")
rm(list=ls())
data <- read_excel("Puertos_data.xlsx")
data <- pdata.frame(data,  index = "port")

yName = "total_traffic"
xNames = c("total_meters", "tugs")
zNames = c("averagehs", "systems_hs")
zIntercept = TRUE
it = c("port", "time")


data(usmanuf)
usmanuf <- pdata.frame(usmanuf,  index = "naics")
get_ml_estimation(yName = "Y",
                  xNames = c("K", "L"),
                  zNames = c("M"),
                  zIntercept = TRUE,
                  it = c("naics", "time"),
                  data = usmanuf)
yName = "Y"
xNames = c("K", "L")
zNames = c("M")
zIntercept = TRUE
it = c("naics", "time")
data = usmanuf

get_ml_estimation = function(yName, xNames, zNames, zIntercept, data, it){
  Y_ = log(data[,yName])
  X_ = data[,xNames]
  Z_ = data[,zNames]
  
  nb = length(X_)
  ng = length(Z_)
  
  for (i in 1:nb){
    assign(paste("beta",i,sep=""), 1)
  }
  if (zIntercept){
    Z_ = cbind.data.frame(replicate(nrow(Z_),0) ,Z_)
    ng = ng + 1
  }
  for (i in 1:ng){
     assign(paste("gamma", i, sep=""), 1)
  }
  
  sigmaesq = 1
  coefficients = c(mget(ls(pattern="beta")), 
                   mget(ls(pattern="gamma")), 
                   mget(ls(pattern="sigmaesq")))
  
  index = unique(data[,"port"])
  time = unique(data[,"time"])
  dataY= cbind.data.frame(Y_, data$port, data$time)
  dataX= cbind.data.frame(X_, data$port, data$time)
  dataZ= cbind.data.frame(Z_, data$port, data$time)
  
  log_likelihood = function(beta1, beta2, gamma1, gamma2, sigmaesq){
    t=length(time)
    get_p_matrix = function(t){
      p = -diag(t)
      p = p[-nrow(p),]
      aux = diag(t)
      aux = aux[-1,]
      p = p+aux
      return(p)
    }
    get_beta_array = function(coefficients, nb){
      beta = array(dim = nb)
      for (i in 1:nb){
        beta[i] = coefficients[[i]]
      }
      return(beta)
    }
    get_gamma_coefficients = function(coefficients, nb, ng){
      gamma = array(dim = ng)
      count = 0
      for (i in (nb+1):(length(coefficients)-1)){
        count = count + 1
        gamma[count] = coefficients[[i]]
      }
      return(gamma)
    }
    
    beta_array = get_beta_array(coefficients, nb)
    gamma_array = get_gamma_coefficients(coefficients, nb, ng)
    sigmaesq = coefficients[[length(coefficients)]]
    
    P = get_p_matrix(t)
    result = array(dim = length(index))
    count = 0
    for (i in index){
      print(i)
      Y = as.matrix(dataY[which(dataY$`data$port`==i),1])
      X = as.matrix(dataX[which(dataX$`data$port`==i),1:nb])
      Z = as.matrix(dataZ[which(dataZ$`data$port`==i),1:ng])
      y = P%*%Y
      x = P%*%X
      R = y - x%*%beta_array
      
      var_array = array(dim=t)
      for (i in 1:t){
        var_array[i]=exp(gamma_array%*%Z[i,])
      }
      
      
      mu = replicate(t-1, 0)
      sigma = as.matrix(P%*%(sigmaesq*diag(t) + diag(var_array))%*%t(P))
      rownames(sigma) <- colnames(sigma)
      gamma = -diag(var_array)%*%t(P)%*%solve(sigma)
      nu = replicate(t, 0)
      delta = diag(var_array)-diag(var_array)%*%t(P)%*%solve(sigma)%*%P%*%diag(var_array)
      
      count = count + 1
      result[count] = dcsn(x=array(R), mu, sigma, gamma, nu, delta)
    }
    print(-sum(log(result)))
    -sum(log(result))
  }
  
  ml = mle2(log_likelihood, start = coefficients,
       method = "L-BFGS-B",
       trace = TRUE,
       lower = c(beta1 =-Inf, beta2=-Inf, gamma1=-Inf, gamma1=-Inf, sigmaesq = 0.001),
       upper = c(beta1 =Inf, beta2=Inf, gamma1=Inf, gamma1=Inf, sigmaesq = 10))
  return(ml)
}






