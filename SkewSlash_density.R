densSSL<-function(x,mu,sigma2,lambda,nu,N){
  y<- matrix(NA,N,1)
  for(j in 1:N){
    f <- function(u){2*nu*u^(nu-1)*dnorm(x[j],mu[j],sqrt(sigma2/u))*pnorm(sqrt(u)*lambda*sqrt(1/sigma2)*(x[j]-mu[j]))}
    y[j,1] <- integrate(f,0,1)$value
  }
  return(y)
}