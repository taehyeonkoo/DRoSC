# setting <- 1
source("~/Documents/GitHub/DRoSC/src/helpers.R")
library(ggplot2)
# library(Matrix)
N <- 10
SC.true <- c(1/3,1/3,1/3,rep(0,N-3))
SC.post <- SC.true
dt <- data.frame()
dt.CI <- data.frame()
# setting <- 2
alpha = 0.05


setting <- as.numeric(readLines(gi,n=1)) # {1,2,3}
T0 <- as.numeric(readLines(gi,n=1)) # {25,50}
T1 <- as.numeric(readLines(gi,n=1))  # {25,50}
phi <- as.numeric(readLines(gi,n=1))  # {0,0.5}
dependent = ifelse(phi>0,T,F)
if (setting==1) { # favorable to SC
  corr <- 0.25
  Sigma0 <- (1-corr)*diag(1,N)+corr*matrix(1,N,N)
  Sigma1 <- diag(1,N)
  mu0 <-   rep(c(0.8,1.2),N/2)#+rnorm(N,0,0.05^2)
  mu1 <- mu0
  # lambda <- 0
  SC.post <- SC.true
} else if(setting ==2){ # X shift, high-correlated
  c <- 0.05
  corr <- 0.95
  Sigma0 <- ((1-corr)*diag(1,N)+corr*matrix(1,N,N)) # bdiag
  
  Sigma1 <- diag(1,N)#Sigma0
  temp <- rnorm(N,0,0.05^2)
  mu0 <-  rep(c(0.8,1.2),N/2)#+rnorm(N,mean = 0, sd = 0.05) 
  
  mu1 <- mu0+c(0.6,0.4,0.2,rep(-0.0,3),rep(-0.0,N-6)) # X shift
  SC.post <- SC.true+c(-rep(c,1),rep(0,N-2),rep(c,1))
  if (any(SC.post<0)|round(sum(SC.post),3)!=1) {
    cat("error for beta1\n")
    break
  }
} else if(setting ==3){ # beta shift
  c <- 0.2
  corr <- 0.25
  scale <- 1
  Sigma0 <- scale^2*((1-corr)*diag(1,N)+corr*matrix(1,N,N))
  Sigma1 <- scale^2*(diag(1,N))
  temp <- rnorm(N,0,0.05^2)
  mu0 <-  1+seq(1,N)/(N)
  mu1 <- mu0
  # c <- as.numeric(readLines(gi,n=1))
  SC.post <- SC.true+c(-rep(c,3),rep(0,N-6),rep(c,3))
  if (any(SC.post<0)|round(sum(SC.post),3)!=1) {
    cat("error for beta1\n")
    break
  }
} 

# T0 <- t0
# T1 <- t1

nsims <- 500
Sigma.true <- Sigma0+mu0%*%t(mu0)
norm.true = norm(Sigma.true,"I") 
# max(abs(Sigma.true%*%(SC.post-SC.true)))
beta.mat <- matrix(nrow = nsims,ncol = N)
mu.mat <- beta.mat
l2.beta <- matrix(nrow = nsims,ncol = N)
linf.beta <- matrix(nrow = nsims,ncol = N)
mu.beta <- matrix(nrow = nsims,ncol = N)
l2.vec <- linf.vec <- mu.vec <- vector()
tau <- as.numeric(readLines(gi,n=1)) # {-1.5,-1.4,...,1.5}
lambda <- max(abs(Sigma.true%*%(SC.post-SC.true)))

if (lambda ==0) {
  beta.star <- SC.post
} else{
  
  g <- rbind(diag(x=1,N,N),Sigma.true, -Sigma.true)
  
  h <- rbind(matrix(0,N,1),cbind(c(Sigma.true%*%SC.true-lambda,
                                   -(Sigma.true%*%SC.true+lambda))))
  beta.star <- limSolve::lsei(A=mu1,B=tau+mu1%*%SC.post,E=matrix(1,1,N),F=1,
                              G=g,H=h,type=2)$X
  
}
tau.star <-c(tau+mu1%*%(SC.post-beta.star))
for (m in 1:nsims) {
  if (m %% 20 ==1) {
    cat("Setting = ",setting,", T0 = ",T0,", T1 = ",T1,", tau = ",tau,", nsim = ",m,"\n")
  }
  F1 <- rnorm(T0+T1)
  F2 <- rnorm(T0+T1)
  eps_y <- generate_AR1_process(T0+T1, 1, phi, 1, 0)
  X0 <- generate_AR1_process(T0, N, phi, Sigma0, mu0)#MASS::mvrnorm(T0,mu = mu0,Sigma = Sigma0)
  
  X1 <-  generate_AR1_process(T1, N, phi, Sigma1, mu1)#MASS::mvrnorm(T1,mu = mu1,Sigma = Sigma1)
  
  mu1hat <- apply(X1,2,mean)
  Y0 <-  as.vector(X0%*%SC.true+eps_y[1:T0])
  Y1 <-  as.vector(X1%*%SC.post+eps_y[-(1:T0)])
  Y1 <- Y1+tau+rnorm(T1,0,sd = 0.25)
  SC.beta <- sc(Y0,X0)$w.hat
  tau.SC <- mean(Y1-X1%*%SC.beta)
  M <- ifelse(phi ==0,500,1000) # M=500 if iid; otherwise M=1000
  
  
  DRoSC.fit <- DRoSC(Y0,Y1,X0,X1,lambda = lambda, dependent = F,c = 1/2,
                   Inference = T,true_beta = beta.star,true_mu = mu1,M = M)
  
  # Sensitivity for order of \rho $
  # DRoSC.fit2 <- DRoSC(Y0,Y1,X0,X1,lambda = lambda, dependent = F,c = 1,
  #                   Inference = F,true_beta = beta.star,true_mu = mu1,M = M)
  # # 
  # # 
  # DRoSC.fit3 <- DRoSC(Y0,Y1,X0,X1,lambda = lambda, dependent = F,c = 2,
  #                   Inference = F,true_beta = beta.star,true_mu = mu1,M = M)
  
  
  dt <- rbind(dt,data.frame(setting = setting,T0 = T0,T1 = T1,tauHat = DRoSC.fit$tauHat,
                            tau = tau,tau.star = tau.star, CI.lower = DRoSC.fit$CI.tau[1],CI.upper = DRoSC.fit$CI.tau[2]))
}
filename <- paste("setting",setting,"-phi",phi,"-pre",T0,"-post",T1,"-tau",tau,".RData", sep="")
save.image(filename)