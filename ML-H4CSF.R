if(!require(csn)){install.packages("csn");library(csn)} 
if(!require(stats4)){install.packages("stats4");library(stats4)} 
if(!require(bbmle)){install.packages("bbmle");library(bbmle)} 
if(!require(matrixcalc)){install.packages("matrixcalc");library(matrixcalc)} 
if(!require(xlsx)){install.packages("xlsx");library(xlsx)}
if(!require(readxl)){install.packages("readxl");library(readxl)} 
if(!require(EnvStats)){install.packages("EnvStats");library(EnvStats)} 
if(!require(plm)){install.packages("plm");library(plm)} 
if(!require(npsf)){install.packages("npsf");library(npsf)} 
if(!require(tmvtnorm)){install.packages("tmvtnorm");library(tmvtnorm)} 


setwd("~/GitHub/data")
rm(list=ls())
data <- read_excel("Puertos_data.xlsx")
data <- pdata.frame(data,  index = "port")

yName = "total_traffic"
xNames = c("total_meters", "tugs", "total_cranes", "storages_total")
zNames = c("averagehs", "averagewind")
zIntercept = TRUE
it = c("port", "time")



#
#data(usmanuf)
#na = (as.character(unique(usmanuf[which(is.na(usmanuf$Y)),"naics"])))
#for (n in na){
#  usmanuf = usmanuf[usmanuf$naics!=n,]
#}
#
#
#usmanuf <- pdata.frame(usmanuf,  index = "naics")
#res = get_ml_estimation(yName = "Y",
#                  xNames = c("K", "L"),
#                  zNames = c("M"),
#                  zIntercept = TRUE,
#                  it = c("naics", "time"),
#                  data = usmanuf)
#yName = "Y"
#xNames = c("K", "L")
#zNames = c("M")
#zIntercept = TRUE
#it = c("naics", "time")
#data = usmanuf


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
  
  log_likelihood = function(beta1, beta2, beta3, beta4, gamma1, gamma2, gamma3, sigmaesq){
    print(paste(beta1, beta2, beta3,gamma1, gamma2, gamma3, sigmaesq))
    t=length(time)
    get_p_matrix = function(t){
      p = -diag(t)
      p = p[-nrow(p),]
      aux = diag(t)
      aux = aux[-1,]
      p = p+aux
      return(p)
    }
    get_mycoast_p_matrix = function(t){
      p = -diag(t)
      p = p[-c((nrow(p)-11):nrow(p)),]
      aux = diag(t)
      aux = aux[-c(1:12),]
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
    beta_array = as.matrix(c(beta1,beta2,beta3,beta4))
    #gamma_array = get_gamma_coefficients(coefficients, nb, ng)
    gamma_array = as.matrix(c(gamma1,gamma2,gamma3))
    #sigmaesq = coefficients[[length(coefficients)]]
    
    P = get_mycoast_p_matrix(t)

    result = array(dim = length(index))

    for (i in 1:length(index)){
      Y = as.matrix(dataY[which(dataY$`data[, it[1]]`==index[i]),1])
      X = as.matrix(dataX[which(dataX$`data[, it[1]]`==index[i]),1:nb])
      Z = as.matrix(dataZ[which(dataZ$`data[, it[1]]`==index[i]),1:ng])
      y = P%*%Y
      x = P%*%X
      R = y - x%*%beta_array
      
      var_array = as.matrix(array(dim=t))
      for (j in 1:t){
        var_array[j]=exp(t(gamma_array)%*%as.matrix(Z[j,]))
      }
      
      var_ui = as.matrix(diag(array(var_array)))
      var_vi = as.matrix(exp(sigmaesq)*diag(t))
      
      mu = replicate(t-12, 0)
      sigma = as.matrix(P%*%(var_vi + var_ui)%*%t(P))
      rownames(sigma)<-colnames(sigma)
      gamma = -var_ui%*%t(P)%*%solve(sigma)
      nu = replicate(t, 0)
      delta = var_ui-var_ui%*%t(P)%*%solve(sigma)%*%P%*%var_ui
      
      #result[count] = dcsn(x=array(R), mu, sigma, gamma, nu, delta)
      #result[i] = loglcsn(x=array(R), mu, sigma, gamma, nu, delta)
      
      density = dmvnorm(x=matrix(R, ncol = length(R)), mean=array(mu), sigma=sigma, log = FALSE)
      distribution = pmvnorm(lower=rep(-Inf,t), upper = as.numeric(array(gamma%*%R)), mean = nu, sigma = delta)
      result[i] = (2**t)*density*distribution[1]
    }
    print(-sum(log(result)))
    -sum(log(result))
  }

  ml = mle2(log_likelihood, start = list(beta1=1,beta2=1,beta3=1,beta4=1,
                                         gamma1=0.001,gamma2=0.001,gamma3=0.001,
                                         sigmaesq=1),
       method = "L-BFGS-B",
       trace = TRUE,
       lower = c(beta1 =-Inf, beta2=-Inf,beta3=-Inf,beta4=-Inf,
                 gamma1=-Inf, gamma2=-Inf,gamma3=-Inf, 
                 sigmaesq = -Inf),
       upper = c(beta1 =Inf, beta2=Inf,beta3=Inf,beta4=Inf,
                 gamma1=Inf, gamma2=Inf,gamma3=Inf, 
                 sigmaesq = +Inf))
  return(ml)
}


res = get_ml_estimation(yName = "total_traffic",
                        xNames = c("total_meters", "tugs", "total_cranes", "storages_total"),
                        zNames = c("averagehs", "averagewind"),
                        zIntercept = TRUE,
                        it = c("port", "time"),
                        data = data)



