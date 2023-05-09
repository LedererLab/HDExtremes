simu_HR <- function(vario_mat, no.simu=1) {
  
  d <- nrow(vario_mat)
  stopifnot(ncol(vario_mat)==d)
  res <- matrix(0, nrow=no.simu, ncol=d)
  cov.mat <- sapply(1:d, function(i) sapply(1:d, function(j) 
    0.5*(vario_mat[i,1] + vario_mat[1,j] - vario_mat[i,j])))
  cov.mat <- cov.mat + 1e-5
  chol.mat <- chol(cov.mat)
  
  for (k in 1:d) {
    trend <- as.vector(vario_mat[,k]/2)
    poisson <- rexp(no.simu)
    bound <- sapply(1:no.simu, function(i) res[i,k])
    while (any(1/poisson > bound)) {
      ind <- (1/poisson > bound)
      n.ind <- sum(ind)
      idx <- (1:no.simu)[ind]
      proc <- t(chol.mat)%*%matrix(rnorm(d*n.ind), ncol=n.ind)
      proc <- exp(t(proc - trend))
      stopifnot(dim(proc)==c(n.ind, d))
      if (k==1) {
        ind.upd <- rep(TRUE, times=n.ind)
      } else {
        ind.upd <- sapply(1:n.ind, function(i) 
          all(1/poisson[idx[i]] * proc[i,1:(k-1)] <= res[idx[i],1:(k-1)]))
      }
      if (any(ind.upd)) {
        idx.upd <- idx[ind.upd]
        res[idx.upd,] <- pmax(res[idx.upd,], 1/poisson[idx.upd]*proc[ind.upd,])
      }
      poisson[ind] <- poisson[ind] + rexp(n.ind)
      bound <- sapply(1:no.simu, function(i) res[i,k])
    } 
  }
  cat("\n")
  return(res)  
}


simu_HRpareto <- function(vario_mat, no.simu=1) {
  
  d <- nrow(vario_mat)
  stopifnot(ncol(vario_mat)==d)
  res <- NULL
  cov.mat <- sapply(1:d, function(i) sapply(1:d, function(j) 
                     0.5*(vario_mat[i,1] + vario_mat[1,j] - vario_mat[i,j])))
  cov.mat <- cov.mat + 1e-5
  chol_mat <- chol(cov.mat)
  while (length(res) < d*no.simu) {
      N <- no.simu - length(res)/d
      shift.idx <- sample(1:d, N, replace=TRUE)
      gauss.proc <- t(chol_mat) %*% matrix(rnorm(d*N), ncol=N)
      gauss.proc <- t(t(gauss.proc) - gauss.proc[cbind(shift.idx,1:N)])
      gauss.proc <- exp(gauss.proc - vario_mat[,shift.idx]/2)
      stopifnot(all(abs(gauss.proc[cbind(shift.idx,1:N)]-1)<=1e-8))
      acc.prob <- apply(gauss.proc, 2, max)/apply(gauss.proc, 2, sum)
      U <- runif(N)
      if (any(U < acc.prob)) {
        res <- cbind(res, gauss.proc[,U<acc.prob])
      }
  }
  res <- 1/(1-runif(no.simu))*t(res)/apply(res,2,max)
  return(res)
}

emp.var <- function(all.data, thresh) {
  
  d <- ncol(all.data)
  n <- nrow(all.data)
  vario_mat <- matrix(0, nrow=d, ncol=d)
  for (i in 1:d) {
    for (j in 1:d) {
      logtrafo_i <- log(1-rank(all.data[,i])/(n+1))
      logtrafo_j <- log(1-rank(all.data[,j])/(n+1))
      for (m in 1:d) {
        ind <- which(all.data[,m] > thresh)
        vario_mat[i,j] <- vario_mat[i,j] + var(logtrafo_i[ind]-logtrafo_j[ind])
      }
    }  
  }
  return(vario_mat/d)
}

myGamma2Theta <- function(Gamma, double.check=TRUE) {
  
  d <- nrow(Gamma)
  stopifnot(all(abs(Gamma-t(Gamma))<1e-5))
  
  Theta <- matrix(0, nrow=d, ncol=d)
  V <- matrix(0, nrow=d, ncol=d)
  V[1,] <- 1
  Sigma <- 0.5*(Gamma%*%V  + t(V)%*%Gamma - Gamma)
  tmp.mu <- 0.5*solve(Sigma[-1,-1], Gamma[-1,1])
  mu     <- c(1-sum(tmp.mu), tmp.mu)
  Theta[-1,-1] <- solve(Sigma[-1,-1])
  Theta[1,] <- -colSums(Theta)
  Theta[,1] <- -rowSums(Theta)
  
  if (double.check) {
    Theta.check <- matrix(0, nrow=d, ncol=d)
    V <- matrix(0, nrow=d, ncol=d)
    V[d,] <- 1
    Sigma <- 0.5*(Gamma%*%V  + t(V)%*%Gamma - Gamma)
    tmp.mu.check <- 0.5*solve(Sigma[-d,-d], Gamma[-d,d])
    mu.check     <- c(tmp.mu.check, 1-sum(tmp.mu.check))
    Theta.check[-d,-d] <- solve(Sigma[-d,-d])
    Theta.check[d,] <- -colSums(Theta.check)
    Theta.check[,d] <- -rowSums(Theta.check)
    stopifnot(all(abs(Theta-Theta.check)<1e-5))
    stopifnot(all(abs(mu-mu.check)<1e-5))
  }
  Lambda <- upper.tri(Theta)*Theta
  diag(Lambda) <- 0
  return(list(Lambda=Lambda, Theta=Theta, mu=mu))
}

Lambda2Theta <- function(Lambda) {
  Theta <- Lambda + t(Lambda)
  diag(Theta) <- - rowSums(Theta)
  return(Theta)
}



################################
###### coordinate descent ######
################################


## objective function with mu' = - (mu-1)
# 4*sum((mu-1)*y) - 4*sum(Theta%*%u) + sum(((mu-1)^2-2*diag(Theta))*diag(v)) - 2*sum((mu-1)*Theta%*%w) + sum(t(Theta)%*%Z%*%Theta)
# = -4*sum(mu'*y) - 4*sum(Theta%*%u) + sum((mu'^2-2*diag(Theta))*diag(v)) + 2*sum(mu'*Theta%*%w) + sum(t(Theta)%*%Z%*%Theta)


library("Rcpp")

cppFunction('NumericMatrix param_update_conditions(NumericVector mu, 
   NumericMatrix Lambda, NumericMatrix u, NumericMatrix v,
   NumericMatrix w, NumericVector y, NumericVector Z, 
   NumericVector tmp_wTheta, NumericMatrix tmp_ZTheta,
   IntegerMatrix indices, int n, double param_Lambda, double param_mu) {
     int d, j, k, l, m;
     d = Lambda.nrow();
     double Zfactor, curr_param, delta_param;
     NumericMatrix new_param(d+1,d);
     for (j = 0; j < d; j++) {
       new_param(j,_) = Lambda(j,_);
     }
     new_param(d,_) = mu;
     for (m=0; m < d*(d+1)/2; m++) {
         j = indices(m,1);
         k = indices(m,0);
         if (k == -1) {
           double tmp = 2*y(j) - tmp_wTheta(j);
           new_param(d,j) = R::sign(tmp)*R::fmax2(fabs(tmp) - sqrt(n)*param_mu, 0)/v(j,j);
         } else {
           Zfactor =  Z(k+k*d+j*d*d) - Z(k+j*d+j*d*d)
                    + Z(j+j*d+k*d*d) - Z(j+k*d+k*d*d)
                    + Z(j+j*d+j*d*d) - Z(j+k*d+j*d*d)
                    + Z(k+k*d+k*d*d) - Z(k+j*d+k*d*d);
           double tmp = 2*u(k,j) + 2*u(j,k) - 2*u(j,j) - 2*u(k,k)
                      +   v(k,j) +   v(j,k) -   v(j,j) -   v(k,k)
                      - new_param(d,j)*(w(k,j)-w(j,j))
                      - new_param(d,k)*(w(j,k)-w(k,k))
                      - tmp_ZTheta(k,j) - tmp_ZTheta(j,k) 
                      + tmp_ZTheta(j,j) + tmp_ZTheta(k,k) 
                      + Zfactor*new_param(k,j);
           curr_param = new_param(k,j);                
           new_param(k,j) = R::sign(tmp)*R::fmax2(fabs(tmp) - sqrt(n)*param_Lambda, 0)/Zfactor;
           delta_param = new_param(k,j) - curr_param;
           if (fabs(delta_param)>0) {
             tmp_wTheta(j) += (w(k,j) - w(j,j))*delta_param;
             tmp_wTheta(k) += (w(j,k) - w(k,k))*delta_param;
             for (l = 0; l < d; l++) {
               tmp_ZTheta(l,j) += (Z(l+k*d+j*d*d) - Z(l+j*d+j*d*d))*delta_param;
               tmp_ZTheta(l,k) += (Z(l+j*d+k*d*d) - Z(l+k*d+k*d*d))*delta_param;
             }
           }
         }
     }
     return new_param; 
   }')


cppFunction('NumericMatrix Zmatrix(NumericMatrix logres) {
  int i, ind, j, l, m, d, n;
  n = logres.nrow();
  d = logres.ncol();
  NumericMatrix output(d*d,d);
  for (j=0; j<d; j++) {
    for (l=0; l<d; l++) {
      for (m=0; m<d; m++) {
        ind = l+m*d;
        output(ind,j) = 0.0;
        for (i=0; i<n; i++) {
          output(ind,j) += logres(i,j)*logres(i,j)*logres(i,l)*logres(i,m);
        }
      }
    }  
  }
  return output;
  }')

cppFunction('NumericMatrix umatrix(NumericMatrix logres) {
  int i, ind, j, l, m, d, n;
  n = logres.nrow();
  d = logres.ncol();
  NumericMatrix output(d,d);
  for (j=0; j<d; j++) {
    for (l=0; l<d; l++) {
      output(l,j) = 0.0;
      for (i=0; i<n; i++) {
        output(l,j) += 0.5*(logres(i,j)*logres(i,j)+2*logres(i,j))*logres(i,l);
      }
    }  
  }
  return output;
  }')

cppFunction('NumericMatrix wmatrix(NumericMatrix logres) {
  int i, ind, j, l, m, d, n;
  n = logres.nrow();
  d = logres.ncol();
  NumericMatrix output(d,d);
  for (j=0; j<d; j++) {
    for (l=0; l<d; l++) {
      output(l,j) = 0.0;
      for (i=0; i<n; i++) {
        output(l,j) += logres(i,j)*logres(i,j)*logres(i,l);
      }
    }  
  }
  return output;
  }')

coord.descent_conditions <- function(res, start.mu, start.Lambda, param_Lambda, param_mu, 
                                     random.order=FALSE, tol=1e-8) {
  
  n <- nrow(res)
  d <- ncol(res)
  stopifnot(length(start.mu)==d)
  stopifnot(all(dim(start.Lambda)==d))
  no.param <- length(param_Lambda)
  stopifnot(length(param_mu)==no.param)
  stopifnot(all(abs(rev(sort(param_Lambda))-param_Lambda) < 1e-10))
  stopifnot(all(abs(rev(sort(param_mu))-param_mu) < 1e-10))  
  res.Lambda <- vector(mode="list", length=no.param)
  res.mu     <- vector(mode="list", length=no.param)
  
  precalc.time <- system.time({
  #### calculate u,v,w,y,Z ####
  u <- umatrix(log(res))  
  v <- diag(apply((log(res))^2, 2, sum)) #former mean
  w <- wmatrix(log(res))
  y <- 0.5*apply( (log(res))^2 + 2*log(res), 2, sum) #former mean
  Z <- Zmatrix(log(res))
  dim(Z) <- c(d,d,d)
  largeZ <- matrix(0, nrow=d^2, ncol=d^2)
  for (j in 1:d) {
    #u[,j] <- 0.5*apply( ((log(res[,j]))^2 + 2*log(res[,j]))*log(res), 2, sum) #former mean
    #w[,j] <- apply((log(res[,j]))^2*log(res), 2, sum) #former mean
    #for (l1 in 1:d) {
    #  for (l2 in 1:d) {
    #    Z[l1,l2,j] <- sum((log(res[,j]))^2*log(res[,l1])*log(res[,l2])) #former mean
    #  }
    #}
     largeZ[(j-1)*d + 1:d, (j-1)*d + 1:d] <- Z[,,j]
  }
  }) 
  
 calc.time <- system.time( {
  #### apply coordinate descent ####
  indices <- as.matrix(expand.grid(-1:(d-1),0:(d-1)))
  indices <- indices[indices[,1] < indices[,2],]
  old.mu <-  mu <- start.mu
  old.Lambda <- Lambda <- start.Lambda
  tmp_wTheta <- colSums(w*Lambda2Theta(old.Lambda))
  tmp_ZTheta <- matrix(NA, nrow=d, ncol=d)
  for (j in 1:d) tmp_ZTheta[,j] <- Z[,,j] %*% Lambda2Theta(old.Lambda)[,j]
  counter <- rep(0, times=no.param)
  param.ind <- 1
  
  while (param.ind <= no.param) {
    if (random.order) {
      indices <- indices[sample(1:(d*(d+1)/2),d*(d+1)/2),]
    }
    counter[param.ind] <- counter[param.ind] + 1
    new_param <- param_update_conditions(old.mu, old.Lambda, u, v, w, y, Z, 
                                         tmp_wTheta, tmp_ZTheta, indices, n, 
                                         param_Lambda[param.ind],
                                         param_mu[param.ind])
    mu <- new_param[d+1,]
    Lambda <- new_param[1:d,]
    counter <- counter + 1
    if (  all(abs(old.Lambda-Lambda) <= tol*pmax(abs(old.Lambda),abs(Lambda))) 
          & all(abs(old.mu-mu)       <= tol*pmax(abs(old.mu),abs(mu))) ) {
      res.Lambda[[param.ind]] <- Lambda
      res.mu[[param.ind]] <- mu
      tmp_wTheta <- colSums(w*Lambda2Theta(Lambda))
      tmp_ZTheta <- matrix(NA, nrow=d, ncol=d)
      for (j in 1:d) tmp_ZTheta[,j] <- Z[,,j] %*% Lambda2Theta(Lambda)[,j]
      #cat(counter, "\n")
      param.ind <- param.ind + 1
    }
    old.mu <- mu
    old.Lambda <- Lambda
  }
 })
  res.Theta <- lapply(res.Lambda, Lambda2Theta)
  stopifnot(length(res.Theta)==length(res.Lambda))
  if (length(res.Lambda)==1) {
    res.Lambda <- res.Lambda[[1]]
    res.Theta  <- res.Theta[[1]]
  }
  return(list(Lambda=res.Lambda, mu=res.mu, Theta=res.Theta, 
              precalc.time=precalc.time[3], calc.time=calc.time[3]))
}