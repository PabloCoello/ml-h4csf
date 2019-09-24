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
na = (as.character(unique(usmanuf[which(is.na(usmanuf$Y)),"naics"])))
for (n in na){
  usmanuf = usmanuf[usmanuf$naics!=n,]
}


usmanuf <- pdata.frame(usmanuf,  index = "naics")
res = get_ml_estimation(yName = "Y",
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
  Y_ = data.frame(log(data[,yName]))
  X_ = data.frame(log(data[,xNames]))
  Z_ = data.frame(data[,zNames])
  
  nb = length(X_)
  ng = length(Z_)
  
  for (i in 1:nb){
    assign(paste("beta",i,sep=""), 1)
  }
  if (zIntercept){
    Z_ = cbind.data.frame(replicate(nrow(Z_),1) ,Z_)
    ng = ng + 1
  }
  for (i in 1:ng){
     assign(paste("gamma", i, sep=""), 0.001)
  }
  
  sigmaesq = 1
  coefficients = c(mget(ls(pattern="beta")), 
                   mget(ls(pattern="gamma")), 
                   mget(ls(pattern="sigmaesq")))
  
  index = unique(data[,it[1]])
  time = unique(data[,"time"])
  dataY= cbind.data.frame(Y_, data[,it[1]], data$time)
  dataX= cbind.data.frame(X_, data[,it[1]], data$time)
  dataZ= cbind.data.frame(Z_, data[,it[1]], data$time)
  
  log_likelihood = function(beta1, beta2, beta3, beta4, beta5, gamma1, gamma2, gamma3, sigmaesq){
    t=length(time)
    get_p_matrix = function(t){
      p = -diag(t)
      p = p[-nrow(p),]
      aux = diag(t)
      aux = aux[-1,]
      p = p+aux
      return(p)
    }
#    get_beta_array = function(coefficients, nb){
#      beta = array(dim = nb)
#      for (i in 1:nb){
#        beta[i] = coefficients[[i]]
#      }
#      return(beta)
#    }
#    get_gamma_coefficients = function(coefficients, nb, ng){
#      gamma = array(dim = ng)
#      count = 0
#      for (i in (nb+1):(length(coefficients)-1)){
#        count = count + 1
#        gamma[count] = coefficients[[i]]
#      }
#      return(gamma)
#    }
    
    #beta_array = get_beta_array(coefficients, nb)
    beta_array = c(beta1,beta2,beta3,beta4,beta5)
    #gamma_array = get_gamma_coefficients(coefficients, nb, ng)
    gamma_array = c(gamma1,gamma2,gamma3)
    #sigmaesq = coefficients[[length(coefficients)]]
    
    P = get_p_matrix(t)
    result = array(dim = length(index))
    count = 0
    for (i in index){
      Y = as.matrix(dataY[which(dataY$`data[, it[1]]`==i),1])
      X = as.matrix(dataX[which(dataX$`data[, it[1]]`==i),1:nb])
      Z = as.matrix(dataZ[which(dataZ$`data[, it[1]]`==i),1:ng])
      y = P%*%Y
      x = P%*%X
      R = y - x%*%beta_array
      
      var_array = array(dim=t)
      for (j in 1:t){
        var_array[j]=exp(gamma_array%*%Z[j,])
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
  
  ml = mle2(log_likelihood, start = list(beta1=1,beta2=1,beta3=1,beta4=1,beta5=1,gamma1=0.001,gamma2=0.001,gamma3=0.001,gamma4=0.001,gamma5=0.001,gamma6=0.001,sigmaesq=1),
       method = "L-BFGS-B",
       trace = TRUE,
       lower = c(beta1 =-Inf, beta2=-Inf,beta3=-Inf,beta4=-Inf,beta5=-Inf, gamma1=-Inf, gamma2=-Inf,gamma3=-Inf, sigmaesq = 0.001),
       upper = c(beta1 =Inf, beta2=Inf,beta3=Inf,beta4=Inf,beta5=Inf, gamma1=Inf, gamma2=Inf,gamma3=Inf, sigmaesq = 10))
  return(ml)
}


res = get_ml_estimation(yName = "total_traffic",
                        xNames = c("total_meters", "tugs", "total_cranes", "storages_total", "wages"),
                        zNames = c("averagehs", "averagewind"),
                        zIntercept = TRUE,
                        it = c("port", "time"),
                        data = data)



