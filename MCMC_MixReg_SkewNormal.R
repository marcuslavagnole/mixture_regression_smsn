#' Bayesian finite mixture of regression models - Skew-Normal
#'
#' This function estimates a finite mixture of regression models, 
#' where each component follows a skew-normal distribution.
#'
#' @param y n-dimensional vector of responses.
#' @param x Matrix - nx(p+1) - of predictors (include the intercept).
#' @param numcomp Number of components cluster to be considered in the 
#' mixture model
#' @param n_mcmc Number of iterations.
#' @param burnin_mcmc Number of initial iterations to be discarded.
#' @param thin_mcmc Thinning parameter.
#' 
#' @return A list with the chains of all parameters of interest.

## Packages
require(FMsmsnReg); require(MCMCpack); require(msm); require(sn) 

## MCMC
mixreg_skewnormal <- function(y,x,numcomp,n_mcmc,burnin_mcmc,thin_mcmc){
  yx<-cbind(y,x)
  yx<-yx[order(yx[,1]),]
  y<-yx[,1]
  x<-yx[,2:(1+dim(x)[2])]
  n <- length(y)
  numcov <- ncol(x)
  resultado <- list()
  # Create auxiliary objects
  mu    <- matrix(NA, n_mcmc, numcomp)
  beta  <- matrix(NA, n_mcmc, numcomp*numcov)
  psi   <- matrix(NA, n_mcmc, numcomp)
  tau   <- matrix(NA, n_mcmc, numcomp)
  p     <- matrix(NA, n_mcmc, numcomp)
  sigma2<- matrix(NA, n_mcmc, numcomp)
  lambda<- matrix(NA, n_mcmc, numcomp)
  s     <- matrix(NA, n_mcmc, n)
  p.aloc<- array(NA, dim = c(numcomp, n_mcmc, n))
  # Set hyperparameter values (hierarchical prior)
  r  <- 1.5                          
  phi<- 0.5                        
  g1 <- 0.5+(r-1)/2                 
  g2 <- g1/(phi*var(y))
  # Set the initial values
  reg_init = regmixEM(y, x[,-1], lambda = NULL, beta = NULL, sigma = NULL, 
                      k = numcomp, addintercept = TRUE, arbmean = TRUE, 
                      arbvar = TRUE, epsilon = 1e-08, maxit = 10000, verb = FALSE)
  beta[1,] <- c(reg_init$beta)
  lambda.i <- rep(0,numcomp) ###
  sigma2.i <- reg_init$sigma^2
  p[1,]    <- reg_init$lambda
  mu[1,]   <- -sqrt(sigma2.i)*(lambda.i/sqrt(1+lambda.i^2))*sqrt(2/pi)
  psi[1,]  <- lambda.i*sqrt(sigma2.i/(1+lambda.i^2))
  tau[1,]  <- sqrt(sigma2.i/(1+lambda.i^2))
  s[1,]    <- sample(1:numcomp,size=n,replace=TRUE,prob=p[1,])
  C        <- rgamma(1,g1,g2)
  #MCMC
  for(k in 2:n_mcmc){
    p[k,] <- atualizarP(s[k-1,],numcomp)
    aux   <- cbind(s[k-1,],y,x)
    for (i in 1:numcomp){
      aux_comp <- subset(aux, aux[,1] == i)
      y_comp   <- aux_comp[,2]
      x_comp   <- aux_comp[,3:(2+numcov)]
      z_comp   <- atualizarZ(y_comp,x_comp,beta[k-1,((i*numcov)-numcov+1):(i*numcov)],mu[k-1,i],tau[k-1,i],psi[k-1,i],nrow(aux_comp))
      x_star   <- cbind(x_comp,z_comp)
      beta.aux <- c(beta[k-1,((i*numcov)-numcov+1):(i*numcov)],psi[k-1,i])
      tau[k,i] <- atualizarTAU(2.5,C,rep(0,numcov+1),diag(c(rep(100,numcov),100)),x_star,beta.aux,mu[k-1,i],y_comp,nrow(aux_comp))
      beta.atualiza <- atualizarBETA(rep(0,numcov+1),diag(c(rep(100,numcov),100)),x_star,mu[k-1,i],tau[k,i],y_comp,nrow(aux_comp))
      beta[k,((i*numcov)-numcov+1):(i*numcov)] <- beta.atualiza[1:numcov] ; psi[k,i] <- beta.atualiza[numcov+1]
      #Original parameters
      sigma2[k,i] <- tau[k,i]^2+psi[k,i]^2
      lambda[k,i] <- psi[k,i]/tau[k,i]
      mu[k,i]     <- -sqrt(sigma2[k,i])*(lambda[k,i]/sqrt(1+lambda[k,i]^2))*sqrt(2/pi)
    }
    C <- atualizarC(g1,g2,2.5,numcomp,tau[k,])
    s.aux <- atualizarS(y,p[k,],mu[k,],x,beta[k,],sigma2[k,],lambda[k,],n,numcomp,numcov)
    s[k,] <- s.aux[,1] ; p.aloc[1:numcomp,k,] = c(s.aux[,-1])
  }
  resultado[['beta']]   <- beta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['sigma2']] <- sigma2[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['lambda']] <- lambda[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['weights']]<- p[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  
  return(resultado)
}


#######################################################################
## Auxiliary functions: Sampling from full conditional distributions ##
#######################################################################

# Allocation vector
atualizarS<-function(y,p,mu,x,beta,sigma2,lambda,N,numcomp,numcov){
  matrizaux<-matrix(NA,numcomp,N)
  for(g in 1:numcomp){
    escala<- x%*%beta[((g*numcov)-numcov+1):(g*numcov)]+mu[g]
    matrizaux[g,]<-p[g]*dsn(y,escala,sqrt(sigma2[g]),lambda[g])
  }
  quociente<-colSums(matrizaux)
  m<-t(matrizaux)/quociente
  cum<-t(apply(m,1,cumsum))
  zero<-matrix(0,N,1)
  acum<-cbind(zero,cum)
  un<-runif(N)
  aux<-ifelse(un>acum,1,0)
  index<-rowSums(aux)
  return(cbind(index,m))
}

# Full conditional for the weights
atualizarP<-function(s,numcomp){
  priori<-matrix(NA,1,numcomp)
  verossi<-matrix(NA,1,numcomp)
  for(j in 1:numcomp){
    priori[1,j]<-4
    verossi[1,j]<-sum(s == j)
  }
  alpha<-priori+verossi
  posteriori<-rdirichlet(1,alpha)
  return(posteriori)
}

# Full conditional for beta
atualizarBETA<-function(b,B,x,mu,tau,y,N){
  B.inv<- chol2inv(chol(B))
  sigma<- chol2inv(chol((1/tau^2)*B.inv+(1/tau^2)*(t(x)%*%x)))
  media<- sigma%*%((1/tau^2)*(B.inv%*%b)+(1/tau^2)*(t(x)%*%(y-mu)))
  beta<- rmvnorm(1,media,sigma)
  return(beta)
}

# Full conditional for tau
atualizarTAU<-function(c,C,b,B,x,beta,mu,y,N){
  alpha <- c+N/2+1/2
  beta <- C+0.5*(t(y-mu-x%*%beta)%*%(y-mu-x%*%beta)+t(beta-b)%*%chol2inv(chol(B))%*%(beta-b))
  tau<-sqrt(1/rgamma(1, alpha, beta))
  return(tau)
}

# Full conditional for the latent variable in asymmetry
atualizarZ<-function(y,x,beta,mu,tau,psi,N){
  media<-(y-x%*%beta-mu)*psi/(tau^2+psi^2)
  sd<-sqrt(tau^2/(tau^2+psi^2))
  z<-rtnorm(N,media,sd,lower=0,upper=Inf)
  return(z)
}

# Full conditional for the hierarchical prior
atualizarC<-function(g1,g2,alpha,numcomp,vetor.tau){
  a<-g1+numcomp*alpha
  b<-g2+sum(1/vetor.tau^2)
  c<-rgamma(1, a, b)
  return(c)
}