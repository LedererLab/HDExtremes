library("evd")
source("functions.R")

n.rep <- 50

HRP.exact <- FALSE
d <- 20
no.simu <- 350000

loc <- 0:(d-1)
floc <- rep(1, times=d) / d^0.25 #cos(pi*loc)

cov.mat <- matrix(NA, nrow=d, ncol=d)
for (i in 1:d) {
  for (j in 1:d) {
    cov.mat[i,j] <- floc[i]*floc[j]*min(loc[i], loc[j])   
  }
}
vario_mat <- matrix(NA, nrow=d, ncol=d)
for (i in 1:d) {
  for (j in 1:d) {
    vario_mat[i,j] <- cov.mat[i,i] - 2*cov.mat[i,j] + cov.mat[j,j]    
  }
}

true_param <- myGamma2Theta(vario_mat)
Theta <- true_param$Theta
Theta <- round(Theta, 10)
Lambda <- true_param$Lambda
Lambda <- round(Lambda, 10)
mu    <- true_param$mu

cat("d=", d, " theoretical non-zeros:", sum(Lambda!=0), "\n")

all_param_Lambda_factors <- c(1000, 100, 10, 1, 0.1, 0.01, 0)
no.param <- length(all_param_Lambda_factors)
all_param_mu <- rep(0, times=no.param)

RMSE_Theta     <- matrix(NA, nrow=no.param, ncol=n.rep)
RMSE_Gamma     <- matrix(NA, nrow=no.param, ncol=n.rep)
no.zeros_Theta <- matrix(NA, nrow=no.param, ncol=n.rep)
precomp.time   <- rep(NA, times=n.rep)
comp.time      <- rep(NA, times=n.rep)

for (k in 1:n.rep) {

cat(k, "")
  
if (HRP.exact) {
  res <- simu_HRpareto(vario_mat=vario_mat, no.simu=no.simu)
} else {
  thresh <- qfrechet(0.95)
  res <- simu_HR(vario_mat=vario_mat, no.simu=no.simu)
  exceed <- apply(res, 1, function(x) max(x) > thresh)
  res <- res[exceed,]/thresh
  print(dim(res))
}
N_u <- nrow(res)  
all_param_Lambda <- all_param_Lambda_factors*sqrt(log(d)/N_u)

tol <- 1e-4

#empvario_mat <- emp.var(all.data=res, thresh=1)
start.mu <- rep(0, times=d)   #myGamma2Theta(empvario_mat)$mu
start.Lambda <- matrix(0, nrow=d, ncol=d) #myGamma2Theta(empvario_mat)$Lambda

  cd.res <- coord.descent_conditions(res=res, start.mu=start.mu, start.Lambda=start.Lambda,
                                     param_Lambda=all_param_Lambda, 
                                     param_mu=all_param_mu)
  est.Theta.list <- cd.res$Theta
  RMSE_Theta[,k] <- unlist(lapply(est.Theta.list, 
                                  function(est.Theta) sqrt(sum((est.Theta-Theta)^2)/d^2)))
  no.zeros_Theta[,k] <- unlist(lapply(est.Theta.list,
                                      function(est.Theta) sum(abs(est.Theta) < 1e-8)))
  precomp.time[k] <- cd.res$precalc.time
  comp.time[k] <- cd.res$calc.time
  
  for (i in 1:no.param) {
    est.Theta <- est.Theta.list[[i]]
    Sigma <- array(0, dim=c(d,d,d))
    for (l in 1:d) {
      Sigma[-l,-l,l] <- solve(est.Theta[-l,-l])
    }
    vario <- array(NA, dim=c(d,d,d))
    for (l in 1:d) {
      for (j in 1:d) {
        vario[l,j,] <- Sigma[l,l,] - 2*Sigma[l,j,] + Sigma[j,j,]
      }
    }
    mean.vario <- apply(vario, 1:2, mean)
    RMSE_Gamma[i,k] <- sqrt(sum((mean.vario-vario_mat)^2)/d^2)
    #var.vario  <- apply(vario, 1:2, var)
  }
  
  #cat("lambda:", param_Lambda, "L2 error lambda:",  sqrt(sum((Theta-cd.res$Theta)^2))/d^2, 
  #    "L2 error Gamma:",  sqrt(sum((mean.vario-vario_mat)^2))/d^2,
  #    "non-zeros:", sum(cd.res$Lambda!=0), "\n\n")
  
}
cat("\n")

cat("precomputation time:", round(mean(precomp.time),3),
                       "(", round(sd(precomp.time),3) , ")\n")
cat("computation time:", round(mean(comp.time),3), 
                    "(", round(sd(comp.time),3), ")\n")
cat("no. of zeros:", round(apply(no.zeros_Theta/(d*(d-1)), 1, mean),3), "\n")
cat("no. of zeros: (", round(apply(no.zeros_Theta/(d*(d-1)), 1, sd),3), ")\n")
cat("RMSE_Theta:", round(apply(RMSE_Theta, 1, mean),3), "\n")
cat("RMSE_Theta: (", round(apply(RMSE_Theta, 1, sd),3), ")\n")
cat("RMSE_Gamma:", round(apply(RMSE_Gamma, 1, mean),3),"\n")
cat("RMSE_Gamma: (", round(apply(RMSE_Gamma, 1, sd),3),")\n")

cat("\n")
########################################
############ here comes CVX ############
########################################

#library("CVXR")
#
#muHat    <- Variable(d)
#ThetaHat <- Variable(rows=d, cols=d, symmetric=TRUE)
#
#u <- matrix(NA, nrow=d, ncol=d)
#v <- diag(apply((log(res))^2, 2, mean))
#y <- 0.5*apply( (log(res))^2 + 2*log(res), 2, mean )
#for (j in 1:d) {
#  u[,j] <- 0.5*apply( ((log(res[,j]))^2 + 2*log(res[,j]))*log(res), 2, mean)
#}
#
#objective <- Minimize(- 2*matrix(y,nrow=1)%*%muHat 
#                      - 2*matrix_trace(u %*% ThetaHat)
#                      - sum(diag(ThetaHat)*diag(v))     
#                      + 0.5/no.simu*sum_squares(
#                          matrix(1,nrow=no.simu,ncol=1)%*%t(muHat)*log(res) + 
#                                        log(res)*log(res)%*%ThetaHat) )
#constraints <- list(sum_entries(ThetaHat, axis=1)==0)
#
#problem <- Problem(objective, constraints)
#result <- solve(problem)
