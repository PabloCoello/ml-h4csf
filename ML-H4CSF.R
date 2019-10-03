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
if(!require(rlist)){install.packages("rlist");library(rlist)} 
if(!require(dplyr)){install.packages("dplyr");library(dplyr)} 
library(readspss)


setwd("~/GitHub/data")
rm(list=ls())
data <- read_excel("Puertos_data.xlsx")
data <- pdata.frame(data,  index = "port")


get_difference_p_matrix = function(t){
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

get_within_p_matrix = function(t){
  I = as.matrix(diag(t))
  l = as.matrix(replicate(t, 1))
  Q = I-1/t*l%*%t(l)
  R = cbind(as.matrix(replicate(t-1, 0)),
            diag(t-1))
  p = R%*%Q
  return(p)
}

get_ll_result = function(mu, sigma, gamma, nu, delta, t, R){
  density = dmvnorm(x=matrix(R, ncol = length(R)), mean=array(mu), sigma=sigma, log = FALSE)
  distribution = pmvnorm(lower=rep(-Inf,t), upper = as.numeric(array(gamma%*%R)), mean = nu, sigma = delta)
  result = (2**t)*density*distribution[1]
  return(result)
}

get_var_array = function(t, gamma_array, Z){
  var_array=exp(as.matrix(Z)%*%gamma_array)
  return(var_array)
}

get_ml_estimation = function(yName, xNames, zNames, zIntercept, data, it, pmatrix){
  Y_ = data.frame(log(data[,yName]))
  X_ = data.frame(log(data[,xNames]))
  Z_ = data.frame(data[,zNames])
  
  nb = length(X_)
  ng = length(Z_)
  
  if (zIntercept){
    Z_ = cbind.data.frame(replicate(nrow(Z_),1) ,Z_)
    ng = ng + 1
  }
  
  index = unique(data[,it[1]])
  time = unique(data[,"time"])
  
  dataY= cbind.data.frame(Y_, data[,it[1]], data$time)
  dataX= cbind.data.frame(X_, data[,it[1]], data$time)
  dataZ= cbind.data.frame(Z_, data[,it[1]], data$time)
  
  log_likelihood = function(beta1, beta2, beta3, beta4, beta5,gamma1, gamma2, gamma3, sigmaesq){
    print(writeLines(c(paste("beta:",beta1, beta2, beta3, beta4, beta5),
                       paste("gamma:", gamma1, gamma2, gamma3),
                       paste("var:", sigmaesq))))
    
    t=length(time)
    beta_array = as.matrix(c(beta1,beta2, beta3, beta4,beta5))
    gamma_array = as.matrix(c(gamma1,gamma2,gamma3))
    
    if (pmatrix=="mycoast"){
      P = get_mycoast_p_matrix(t)
      tn = 12
    }else if(pmatrix=="within"){
      P = get_within_p_matrix(t)
      tn = 1
    }else if(pmatrix=="difference"){
      P = get_difference_p_matrix(t)
      tn = 1
    }else{
      stop("P matrix not defined")
    }
    
    result = array(dim = length(index))
    result1 = array(dim = length(index))
    result2 = array(dim = length(index))
    j=0
    time <- proc.time()
    for (i in 1:length(index)){
      Y = as.matrix(dataY[which(dataY$`data[, it[1]]`==index[i]),1])
      X = as.matrix(dataX[which(dataX$`data[, it[1]]`==index[i]),1:nb])
      Z = as.matrix(dataZ[which(dataZ$`data[, it[1]]`==index[i]),1:ng])
      
      y = P%*%Y
      x = P%*%X
      
      R = y - (x%*%beta_array)
      
      var_array = get_var_array(t, gamma_array, Z)
      
      var_ui = as.matrix(diag(array(var_array)))
      var_vi = as.matrix(sigmaesq*diag(t))
      
      mu = replicate(t-tn, 0)
      sigma = as.matrix(P%*%(var_vi + var_ui)%*%t(P));rownames(sigma)<-colnames(sigma)
      gamma = -(var_ui%*%t(P)%*%solve(sigma))
      nu = replicate(t, 0)
      delta = var_ui-(var_ui%*%t(P)%*%solve(sigma))%*%(P%*%var_ui)
      
      result[i] = get_ll_result(mu, sigma, gamma, nu, delta, t, R)
      #result[i] = dcsn(x=array(R), mu, sigma, gamma, nu, delta)
      #result2[i] = loglcsn(x=array(R), mu, sigma, gamma, nu, delta)
    }
    for (i in 1:length(result)){
      if (result[i] == 0){
        result[i]=1e-200
      }
    }
    j=j+1
    print(proc.time()-time)
    print(paste("m3:|||||",-sum(log(result)),"|||||"))
    print(paste("number of iterations:",(j)))
    print("MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM")
    
    -sum(log(result))
  }
  ols = lm(dataY[,1] ~ dataX[,1]+dataX[,2]+dataX[,3]+dataX[,4]+dataX[,5])
  ml = mle2(log_likelihood, start = list(beta1 = as.numeric(ols$coefficients[2]),
                                         beta2 = as.numeric(ols$coefficients[3]),
                                         beta3 = as.numeric(ols$coefficients[4]),
                                         beta4 = as.numeric(ols$coefficients[5]),
                                         beta5 = as.numeric(ols$coefficients[6]),
                                         gamma1=0.5,gamma2=0.5,gamma3=0.5,
                                         sigmaesq=0.0217331999366858),
            method = "L-BFGS-B",
            lower = c(beta1 =-Inf, beta2=-Inf,beta3 =-Inf, beta4=-Inf,beta5=-Inf,
                      gamma1=-Inf, gamma2=-Inf,gamma3=-Inf,
                      sigmaesq = 0.0001),
            upper = c(beta1 =Inf, beta2=Inf,beta3 =Inf, beta4=Inf,beta5=Inf,
                      gamma1=Inf, gamma2=Inf,gamma3=Inf,
                      sigmaesq = 10))
  
  return(ml)
}

res = get_ml_estimation(yName = "total_traffic",
                        xNames = c("total_meters", "tugs", "total_cranes", "wages", "storages_total"),
                        zNames = c("averagehs", "averagewind"),
                        zIntercept = TRUE,
                        it = c("port", "time"),
                        data = data,
                        pmatrix = "within" )#c("mycoast", "within", "difference")



install.packages("drat")
drat::addRepo("JanMarvin")
install.packages("readspss")

library(readspss)

log_likelihood_nls = function ()
var_nui = exp(W %*% delta_array)
var_uit = Z %*% gamma_array

er_it = beta0 -sqrt(2/pi)*(sqrt(var_nui)+ sqrt(var_uit))


er2_it = var_tau + var_nu + var_v + var_u + (beta0 - sqrt(2/pi)*(sqrt(var_nu)+ sqrt(var_u))**2)
err_itis


model = list(
  er = beta0 -sqrt(2/pi)*(sqrt(var_nui)+ sqrt(var_uit))
  er2 = var_tau + var_nu + var_v + var_u + (beta0 - sqrt(2/pi)*(sqrt(var_nu)+ sqrt(var_u))**2)
  err
)

nlsur(eqns=model, data=, type="FGNLS")
exp(10)

