#' Bayesian finite mixture of regression models - Skew-Slash
#'
#' This function estimates a finite mixture of regression models, 
#' where each component follows a skew-slash distribution.
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
mixreg_skewslash <- function(y,x,numcomp,n_mcmc,burnin_mcmc,thin_mcmc){
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
  nu    <- matrix(NA, n_mcmc, numcomp)
  sigma2<- matrix(NA, n_mcmc, numcomp)
  lambda<- matrix(NA, n_mcmc, numcomp)
  s     <- matrix(NA, n_mcmc, n)
  p.aloc<- array(NA, dim = c(numcomp, n_mcmc, n))
  # Create auxiliary objects for the adaptative MH
  rmw <- matrix(NA, n_mcmc, 3*numcomp)
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
  nu[1,]   <- rep(3,numcomp) ###
  mu[1,]   <- -sqrt(sigma2.i)*(lambda.i/sqrt(1+lambda.i^2))*sqrt(2/pi)
  psi[1,]  <- lambda.i*sqrt(sigma2.i/(1+lambda.i^2))
  tau[1,]  <- sqrt(sigma2.i/(1+lambda.i^2))
  s[1,]    <- sample(1:numcomp,size=n,replace=TRUE,prob=p[1,])
  C        <- rgamma(1,g1,g2)
  # Set the initial values for the adaptative MH
  rmw[1,]  <- rep(c(0.8,1,0),numcomp)
  #
  aux     <- cbind(s[1,],y)
  dataaux <- NULL
  uaux    <-NULL
  for (i in 1:numcomp){
    x_comp  <- subset(aux, aux[,1] == i)
    dataaux <- c(dataaux, x_comp[,2])
    uaux    <- c(uaux,rbeta(length(x_comp[,1]),nu[1,i],1))
  }
  v1 <- data.frame(x1=y)
  v2 <- data.frame(x1=dataaux, x2=uaux)
  v  <- unique(merge(v1,v2,by="x1"))
  #MCMC
  for(k in 2:n_mcmc){
    p[k,]    <- atualizarP(s[k-1,],numcomp)
    aux      <- cbind(s[k-1,],v,x)
    dataaux <- NULL
    uaux     <- NULL
    for (i in 1:numcomp){
      aux_comp <- as.matrix(subset(aux, aux[,1] == i))
      y_comp   <- aux_comp[,2]
      x_comp   <- aux_comp[,4:(3+numcov)]
      z_comp   <- atualizarZ(y_comp,x_comp,beta[k-1,((i*numcov)-numcov+1):(i*numcov)],mu[k-1,i],tau[k-1,i],psi[k-1,i],aux_comp[,3],nrow(aux_comp))
      u_comp   <- atualizarU(nu[k-1,i],y_comp,x_comp,beta[k-1,((i*numcov)-numcov+1):(i*numcov)],mu[k-1,i],psi[k-1,i],tau[k-1,i],z_comp,nrow(aux_comp))
      x_star   <- cbind(x_comp,z_comp)
      beta.aux <- c(beta[k-1,((i*numcov)-numcov+1):(i*numcov)],psi[k-1,i])
      tau[k,i] <- atualizarTAU(2.5,C,rep(0,numcov+1),diag(c(rep(100,numcov),100)),x_star,beta.aux,mu[k-1,i],u_comp,y_comp,nrow(aux_comp))
      beta.atualiza <- atualizarBETA(rep(0,numcov+1),diag(c(rep(100,numcov),100)),x_star,mu[k-1,i],tau[k,i],u_comp,y_comp,nrow(aux_comp))
      beta[k,((i*numcov)-numcov+1):(i*numcov)] <- beta.atualiza[1:numcov] ; psi[k,i]<-beta.atualiza[numcov+1]
      nu.aux   <- atualizarNU(nu[k-1,i],u_comp,nrow(aux_comp),rmw[k-1,(1+3*(i-1))],rmw[k-1,(2+3*(i-1))],rmw[k-1,(3+3*(i-1))],k)
      nu[k,i] <- nu.aux[1] ; rmw[k,(1+3*(i-1)):(3+3*(i-1))] <- nu.aux[3:5]
      dataaux  <- c(dataaux,y_comp)
      uaux     <- c(uaux,u_comp)
      # Original parameters
      sigma2[k,i] <- tau[k,i]^2+psi[k,i]^2
      lambda[k,i] <- psi[k,i]/tau[k,i]
      mu[k,i]     <- -sqrt(sigma2[k,i])*(lambda[k,i]/sqrt(1+lambda[k,i]^2))*sqrt(nu[k,i]/pi)*(gamma((nu[k,i]-1)/2)/gamma(nu[k,i]/2))
    }
    C <- atualizarC(g1,g2,2.5,numcomp,tau[k,])
    v1<-data.frame(x1=y)
    v2<-data.frame(x1=dataaux,x2=uaux)
    v <-unique(merge(v1,v2,by="x1"))
    s.aux <- atualizarS(y,p[k,],mu[k,],x,beta[k,],sigma2[k,],lambda[k,],nu[k,],n,numcomp,numcov)
    s[k,] <- s.aux[,1] ; p.aloc[1:numcomp,k,] = c(s.aux[,-1])
  }
  resultado[['beta']]   <- beta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['sigma2']] <- sigma2[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['lambda']] <- lambda[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['nu']]     <- nu[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['weights']]<- p[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  
  return(resultado)
}


#######################################################################
## Auxiliary functions: Sampling from full conditional distributions ##
#######################################################################

# Allocation vector
atualizarS<-function(y,p,mu,x,beta,sigma2,lambda,nu,N,numcomp,numcov){
  matrizaux<-matrix(NA,numcomp,N)
  for(g in 1:numcomp){
    escala<- x%*%beta[((g*numcov)-numcov+1):(g*numcov)]+mu[g]
    matrizaux[g,]<-p[g]*densSSL(y,escala,sigma2[g],lambda[g],nu[g],N)
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
atualizarBETA<-function(b,B,x,mu,tau,u,y,N){
  B.inv<- chol2inv(chol(B))
  sigma<- chol2inv(chol((1/tau^2)*B.inv+(1/tau^2)*(t(u*x)%*%x)))
  media<- sigma%*%((1/tau^2)*(B.inv%*%b)+(1/tau^2)*(t(u*x)%*%(y-mu)))
  beta<- rmvnorm(1,media,sigma)
  return(beta)
}

# Full conditional for tau
atualizarTAU<-function(c,C,b,B,x,beta,mu,u,y,N){
  alpha <- c+N/2+1/2
  beta <- C+0.5*(t(u*(y-mu-x%*%beta))%*%(y-mu-x%*%beta)+t(beta-b)%*%chol2inv(chol(B))%*%(beta-b))
  tau<-sqrt(1/rgamma(1, alpha, beta))
  return(tau)
}

# Full conditional for the latent variable in asymmetry
atualizarZ<-function(y,x,beta,mu,tau,psi,u,N){
  media<-(y-mu-x%*%beta)*psi/(tau^2+psi^2)
  sd<-sqrt(tau^2/(u*(tau^2+psi^2)))
  z<-rtnorm(N,media,sd,lower=0,upper=Inf)
  return(z)
}

# Full conditional for the latent variable in nu
atualizarU<-function(nu,y,x,beta,mu,psi,tau,z,N){
  alpha<-c(rep(nu+1,N))
  beta<-(y-mu-x%*%beta-psi*z)^2/(2*tau^2)+0.5*z^2
  aux<-runif(N)
  u<-qgamma(aux*pgamma(1,alpha,beta),alpha,beta)
}

# Full conditional for nu 
condicionalNU<-function(nu,u,N){
  priori<-5*log(nu)-(1*nu)
  verossi<- N*log(nu)+nu*sum(log(u))
  funcao<-priori+verossi
  return(funcao)
}

# Metropolis-Hastings for nu
atualizarNU<-function(nu,u,N,clap,clap.aux,M0,t){
  valoratual<-nu
  valorproposto<-rtnorm(1,valoratual,sqrt(clap*clap.aux),lower=1.4,upper=40)
  candidato<-exp(condicionalNU(valorproposto,u,N)-condicionalNU(valoratual,u,N)-dtnorm(valorproposto,valoratual,sqrt(clap*clap.aux),lower=1.4,upper=40,log=TRUE)+dtnorm(valoratual,valorproposto,sqrt(clap*clap.aux),lower=1.4,upper=40,log=TRUE))
  chanceaceitar<-min(1,candidato)
  contador<-NULL
  if (runif(1)<chanceaceitar){
    NUfinal<-valorproposto
    contador<-1
  } 
  else{
    NUfinal<-valoratual 
    contador<-0
  }
  gama1<- 1/t^0.5
  gama2<- 1*gama1
  termometro<- exp(log(clap)+gama2*(chanceaceitar-0.234))
  termometro.aux<- clap.aux+gama1*((NUfinal-M0)^2-clap.aux)
  p.auxiliar<-M0+gama2*(NUfinal-M0)
  return(c(NUfinal,contador,termometro,termometro.aux,p.auxiliar))
}

# Full conditional for the hierarchical prior
atualizarC<-function(g1,g2,alpha,numcomp,vetor.tau){
  a<-g1+numcomp*alpha
  b<-g2+sum(1/vetor.tau^2)
  c<-rgamma(1, a, b)
  return(c)
}


##################################################
## Auxiliary functions: Supplementary functions ##
##################################################

densSSL<-function(x,mu,sigma2,lambda,nu,N){
  y<- matrix(NA,N,1)
  for(j in 1:N){
    f <- function(u){2*nu*u^(nu-1)*dnorm(x[j],mu[j],sqrt(sigma2/u))*pnorm(sqrt(u)*lambda*sqrt(1/sigma2)*(x[j]-mu[j]))}
    y[j,1] <- integrate(f,0,1)$value
  }
  return(y)
}