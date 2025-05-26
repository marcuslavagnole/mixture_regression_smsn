# Sample the allocation vector
atualizarS<-function(dados,p,mu,covariaveis,beta,sigma2,lambda,N,numcomp,numcov){
  matrizaux<-matrix(NA,numcomp,N)
  for(g in 1:numcomp){
    escala<- covariaveis%*%beta[((g*numcov)-numcov+1):(g*numcov)]+mu[g]
    matrizaux[g,]<-p[g]*dsn(dados,escala,sqrt(sigma2[g]),lambda[g])
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
atualizarBETA<-function(b,B,x,mu,tau,dados,N){
  B.inv<- chol2inv(chol(B))
  sigma<- chol2inv(chol((1/tau^2)*B.inv+(1/tau^2)*(t(x)%*%x)))
  media<- sigma%*%((1/tau^2)*(B.inv%*%b)+(1/tau^2)*(t(x)%*%(dados-mu)))
  beta<- rmvnorm(1,media,sigma)
  return(beta)
}

# Full conditional for tau
atualizarTAU<-function(c,C,b,B,x,beta,mu,dados,N){
  alpha <- c+N/2+1/2
  beta <- C+0.5*(t(dados-mu-x%*%beta)%*%(dados-mu-x%*%beta)+t(beta-b)%*%chol2inv(chol(B))%*%(beta-b))
  tau<-sqrt(1/rgamma(1, alpha, beta))
  return(tau)
}

# Full conditional for the latent variable in asymmetry
atualizarZ<-function(dados,x,beta,mu,tau,psi,N){
  media<-(dados-x%*%beta-mu)*psi/(tau^2+psi^2)
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
