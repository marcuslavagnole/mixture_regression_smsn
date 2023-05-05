# Load libraries
library(coda)
library(gtools)
library(msm)
library(mvtnorm)
library(sn)
library(mixtools)

# Set working directory to file location and load utils
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Full_Conditionals/FullConditionals_MixReg_SkewSlash.R")
source("SkewSlash_density.R")
data(tonedata)

# Define dependent and indepent variables
covariaveis<-cbind(rep(1,length(tonedata[,1])),tonedata[,1])
dados<-tonedata[,2]
yx<-cbind(dados,covariaveis)
yx<-yx[order(yx[,1]),]
dados<-yx[,1]
covariaveis<-yx[,2:(1+dim(covariaveis)[2])]
numcov<-dim(covariaveis)[2]

# Set hyperparameters values
b<-rep(0,numcov+1)              # beta
B<-diag(c(rep(100,numcov),100)) # beta
c<- 2.5                         # tau
r<-1.5                          # hierarchical prior
phi<-0.5                        # hierarchical prior
g1<-0.5+(r-1)/2                 # hierarchical prior
g2<-g1/(phi*var(dados))         # hierarchical prior

# Define number of iterations and components  
NN=20000
numcomp<-2

# Create auxiliary objects
mu <- matrix(NA, NN, numcomp)
beta<- matrix(NA, NN, numcomp*numcov)
psi <- matrix(NA, NN, numcomp)
tau<-matrix(NA, NN, numcomp)
p<- matrix(NA, NN, numcomp)
nu<-matrix(NA, NN, numcomp)
sigma2<-matrix(NA, NN, numcomp)
lambda<-matrix(NA, NN, numcomp)
s<- matrix(NA, NN, length(dados))
p1<- matrix(NA, NN, length(dados))
p2<- matrix(NA, NN, length(dados))

u<-list()
z<-list()
y<-list()

# Set the initial values
beta[1,]<-c(1.9103,0.0354,0.0019,0.9978)
lambda.i<- c(-8.0653,-0.7363)
sigma2.i<-c(0.0029,0.0001)
p[1,]<-c(0.5496,(1-0.5496))
nu[1,]<-c(2.1,2.1)

mu[1,]<- -sqrt(sigma2.i)*(lambda.i/sqrt(1+lambda.i^2))*sqrt(2/pi)*(gamma((nu[1,]-1)/2)/gamma(nu[1,]/2))
psi[1,]<-lambda.i*sqrt(sigma2.i/(1+lambda.i^2))
tau[1,]<-sqrt(sigma2.i/(1+lambda.i^2))
s[1,]<-sample(1:numcomp,size=length(dados),replace=TRUE,prob=p[1,])
C<-rgamma(1,g1,g2)

aux<-cbind(s[1,],dados)
dadosaux<-NULL
uaux<-NULL
for (i in 1:numcomp){
  x<-subset(aux, aux[,1] == i)
  dadosaux<-c(dadosaux,x[,2])
  uaux<-c(uaux,rbeta(length(x[,1]),nu[1,i],1))
}
v1<-data.frame(x1=dados)
v2<-data.frame(x1=dadosaux,x2=uaux)
v<-unique(merge(v1,v2,by="x1"))

# Create auxiliary objects for the adaptative MH
rmw<-matrix(NA,NN,6)
rmw[1,]<-c(0.8,1,0,0.8,1,0)
contador<-rep(0,numcomp)

#MCMC
for(k in 2:NN){
  p[k,]<-atualizarP(s[k-1,],numcomp)
  aux<-cbind(s[k-1,],v,covariaveis)
	dadosaux<-NULL
	uaux<-NULL
	for (i in 1:numcomp){
	  y[[i]]<-as.matrix(subset(aux, aux[,1] == i))
		z[[i]]<-atualizarZ(y[[i]][,2],y[[i]][,4:(3+numcov)],beta[k-1,((i*numcov)-numcov+1):(i*numcov)],mu[k-1,i],tau[k-1,i],psi[k-1,i],y[[i]][,3],length(y[[i]][,1]))
		u[[i]]<-atualizarU(nu[k-1,i],y[[i]][,2],y[[i]][,4:(3+numcov)],beta[k-1,((i*numcov)-numcov+1):(i*numcov)],mu[k-1,i],psi[k-1,i],tau[k-1,i],z[[i]],length(y[[i]][,1]))
		x<-cbind(y[[i]][,4:(3+numcov)],z[[i]])
		beta.aux<-c(beta[k-1,((i*numcov)-numcov+1):(i*numcov)],psi[k-1,i])
    tau[k,i]<-atualizarTAU(c,C,b,B,x,beta.aux,mu[k-1,i],u[[i]],y[[i]][,2],length(y[[i]][,1]))
    beta.atualiza<-atualizarBETA(b,B,x,mu[k-1,i],tau[k,i],u[[i]],y[[i]][,2],length(y[,1]));beta[k,((i*numcov)-numcov+1):(i*numcov)]<-beta.atualiza[1:numcov];psi[k,i]<-beta.atualiza[numcov+1]
    nu.aux<-atualizarNU(nu[k-1,i],u[[i]],length(y[[i]][,1]),rmw[k-1,(1+3*(i-1))],rmw[k-1,(2+3*(i-1))],rmw[k-1,(3+3*(i-1))],k);nu[k,i]<-nu.aux[1];contador[i]<-contador[i]+nu.aux[2];rmw[k,(1+3*(i-1)):(3+3*(i-1))]<-nu.aux[3:5]
		dadosaux<-c(dadosaux,y[[i]][,2])
		uaux<-c(uaux,u[[i]])
		#Original parameters
		sigma2[k,i]<- tau[k,i]^2+psi[k,i]^2
		lambda[k,i]<- psi[k,i]/tau[k,i]
		mu[k,i]<- -sqrt(sigma2[k,i])*(lambda[k,i]/sqrt(1+lambda[k,i]^2))*sqrt(2/nu[k,i])*(2*nu[k,i]/(2*nu[k,i]-1))
	}
	C<-atualizarC(g1,g2,c,numcomp,tau[k,])
	v1<-data.frame(x1=dados)
	v2<-data.frame(x1=dadosaux,x2=uaux)
	v<-unique(merge(v1,v2,by="x1"))
	s.aux<-atualizarS(dados,p[k,],mu[k,],covariaveis,beta[k,],sigma2[k,],lambda[k,],nu[k,],length(dados),numcomp,numcov);s[k,]<-s.aux[,1];p1[k,]<-s.aux[,2];p2[k,]<-s.aux[,3]        
}
