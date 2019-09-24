if(!require(csn)){install.packages("csn");library(csn)} 
if(!require(stats4)){install.packages("stats4");library(stats4)} 
if(!require(bbmle)){install.packages("bbmle");library(bbmle)} 
if(!require(matrixcalc)){install.packages("matrixcalc");library(matrixcalc)} 
if(!require(xlsx)){install.packages("xlsx");library(xlsx)}
if(!require(readxl)){install.packages("readxl");library(readxl)} 
if(!require(EnvStats)){install.packages("EnvStats");library(EnvStats)} 
if(!require(plm)){install.packages("plm");library(plm)} 

setwd("~/GitHub/data")
rm(list=ls())
data <- read_excel("Puertos_data.xlsx")
data <- pdata.frame(data,  index = "port")

yName = "total_traffic"
xNames = c("total_meters", "tugs")
zNames = c("averagehs", "systems_hs")
zIntercept = TRUE
it = c("port", "time")
get_ml_estimation(yName = "total_traffic",
                  xNames = c("total_meters", "tugs"),
                  zNames = c("averagehs", "systems_hs"),
                  zIntercept = TRUE,
                  it = c("port", "time"),
                  data = data)


get_ml_estimation = function(yName, xNames, zNames, zIntercept, data, it){
  Y = log(data[,yName])
  X = data[,xNames]
  Z = data[,zNames]
  
  for (i in 1:length(X)){
    assign(paste("beta",i,sep=""), 1)
  }
  if (zIntercept){
    for (i in 1:length(Z)+1){
      assign(paste("gamma", i, sep=""), 1)
    }
    Z = cbind.data.frame(replicate(nrow(Z),0) ,Z)
  }else{
    for (i in 1:length(Z)){
      assign(paste("gamma", i, sep=""), 1)
    }
  }
  sigmaesq = 1
  coefficients = c(mget(ls(pattern="gamma")), 
                   mget(ls(pattern="beta")), 
                   mget(ls(pattern="sigmaesq")))
  
  index = unique(data[,"port"])
  time = unique(data[,"time"])
  dataY= cbind.data.frame(Y, data$port, data$time)
  dataX= cbind.data.frame(X, data$port, data$time)
  dataZ= cbind.data.frame(Z, data$port, data$time)
  
  log_likelihood = function(...){
    t=length(time)
    get_p_matrix = function(t){
      p = -diag(t)
      p = p[-nrow(p),]
      aux = diag(t)
      aux = aux[-1,]
      p = p+aux
      return(p)
    }
    get_beta_array = function(coefficients, X){
      beta = array(dim = length(X))
      for (i in 1:length(X)){
        beta[i] = coefficients[[i]]
      }
      return(beta)
    }
    get_gamma_coefficients = function(coefficients, X, Z){
      gamma = array(dim = length(Z))
      for (i in length(X):length(Z)){
        gamma[i] = coefficients[[i]]
      }
      return(gamma)
    }
    
    beta = get_beta_array(coefficients, X)
    gamma = get_gamma_coefficients(coefficients, X, Z)
    sigmaesq = coefficients[[length(coefficients)]]
    
    P = get_p_matrix(t)
    result = array(dim = length(index))
    count = 0
    for (i in index){
      Y = as.matrix(dataY[which(dataY$`data$port`==i),1])
      X = as.matrix(dataX[which(dataX$`data$port`==i),1:length(beta)])
      Z = as.matrix(dataZ[which(dataZ$`data$port`==i),1:length(gamma)])
      y = P%*%Y
      x = P%*%X
      R = y - x%*%beta
      
      var_array = array(dim=t)
      for (i in 1:t){
        var_array[i]=exp(gamma%*%Z[i])
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

    -sum(log(result))
  }
  
  ml = mle2(log_likelihood(coefficients), start = coefficients,
       method = "L-BFGS-B",
       trace = TRUE)
  return(ml)
}






