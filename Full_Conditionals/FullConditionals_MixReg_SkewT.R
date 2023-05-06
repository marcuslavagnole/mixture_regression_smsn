# Sample the allocation vector
atualizarS<-function(dados,p,mu,covariaveis,beta,sigma2,lambda,nu,N,numcomp,numcov){
  matrizaux<-matrix(NA,numcomp,N)
  for(g in 1:numcomp){
    escala<- covariaveis%*%beta[((g*numcov)-numcov+1):(g*numcov)]+mu[g]
    matrizaux[g,]<-p[g]*dst(dados,escala,sqrt(sigma2[g]),lambda[g],nu=nu[g])
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
atualizarBETA<-function(b,B,x,mu,tau,u,dados,N){
  B.inv<- chol2inv(chol(B))
  sigma<- chol2inv(chol((1/tau^2)*B.inv+(1/tau^2)*(t(u*x)%*%x)))
  media<- sigma%*%((1/tau^2)*(B.inv%*%b)+(1/tau^2)*(t(u*x)%*%(dados-mu)))
  beta<- rmvnorm(1,media,sigma)
  return(beta)
}
# Full conditional for tau
atualizarTAU<-function(c,C,b,B,x,beta,mu,u,dados,N){
  alpha <- c+N/2+1/2
  beta <- C+0.5*(t(u*(dados-mu-x%*%beta))%*%(dados-mu-x%*%beta)+t(beta-b)%*%chol2inv(chol(B))%*%(beta-b))
  tau<-sqrt(1/rgamma(1, alpha, beta))
  return(tau)
}
# Full conditional for latent variable in asymmetry
atualizarZ<-function(dados,x,beta,mu,tau,psi,u,N){
  media<-(dados-mu-x%*%beta)*psi/(tau^2+psi^2)
  sd<-sqrt(tau^2/(u*(tau^2+psi^2)))
  z<-rtnorm(N,media,sd,lower=0,upper=Inf)
  return(z)
}
# Full conditional for the latent variable in the degrees of freedom
atualizarU<-function(nu,dados,x,beta,mu,psi,tau,z,N){
  alpha<-c(rep(nu/2+1,N))
  beta<-nu/2+(dados-mu-x%*%beta-psi*z)^2/(2*tau^2)+0.5*z^2
  u<-rgamma(N,alpha,beta)
  return(u)
}
# Full conditional for the degrees of freedom 
condicionalNU<-function(nu,u,N){
  d<-4/(1+sqrt(2))
  priori<- log(nu)-3*log(nu+d)
  verossi<- (N*nu/2)*log(nu/2)-N*log(gamma(nu/2))+(nu/2-1)*sum(log(u))-nu/2*sum(u)
  funcao<-priori+verossi
  return(funcao)
}
# Metropolis-Hasting for the degrees of freedom 
atualizarNU<-function(nu,u,N,clap,clap.aux,M0,t){
  valoratual<-nu
  valorproposto<-rtnorm(1,valoratual,sqrt(clap*clap.aux),lower=2,upper=40)
  candidato<-exp(condicionalNU(valorproposto,u,N)-condicionalNU(valoratual,u,N)-dtnorm(valorproposto,valoratual,sqrt(clap*clap.aux),lower=2,upper=40,log=TRUE)+dtnorm(valoratual,valorproposto,sqrt(clap*clap.aux),lower=2,upper=40,log=TRUE))
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