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


source("~/ARSC/src/helpers.R")
# library(ggplot2)
# library(Matrix)

gi = file("stdin")
open(gi)

N <- 10
SC.true <- c(1/3,1/3,1/3,rep(0,N-3))
SC.post <- SC.true
dt <- data.frame()
dt.CI <- data.frame()
setting <- 3
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
  mu0 <-   rep(c(0.8,1.2),N/2)+rnorm(N,0,0.05^2)
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
  c <- 1/6
  corr <- 0.25
  scale <- 1
  Sigma0 <- scale^2*((1-corr)*diag(1,N)+corr*matrix(1,N,N))
  Sigma1 <- scale^2*(diag(1,N))
  temp <- rnorm(N,0,0.05^2)
  mu0 <-  seq(1,N)/(N)
  mu1 <- mu0
  # c <- as.numeric(readLines(gi,n=1))
  SC.post <- SC.true+c(-rep(c,3),rep(0,N-6),rep(c,3))
  if (any(SC.post<0)|round(sum(SC.post),3)!=1) {
    cat("error for beta1\n")
    break
  }
} 
sqrt(mu0^2)
0.025*max(eigen(t(X0)%*%X0/T0)$values)
0.1*min(1/maxnorm)
# T0 <- t0
# T1 <- t1

nsims <- 500
Sigma.true <- Sigma0+mu0%*%t(mu0)
norm.true = norm(Sigma.true,"I") 
# max(abs(Sigma.true%*%(SC.post-SC.true)))
beta.mat <- matrix(nrow = nsims,ncol = N)
l2.beta <- matrix(nrow = nsims,ncol = N)
linf.beta <- matrix(nrow = nsims,ncol = N)
mu.beta <- matrix(nrow = nsims,ncol = N)
l2.vec <- linf.vec <- mu.vec <- vector()
# i <- 1
# i <- 1
tau <- as.numeric(readLines(gi,n=1))
# tau <- -2
# star.vec <- vector()
# temp <- apply(X0,2,mean)
lambda <- max(abs(Sigma.true%*%(SC.post-SC.true)))#1.75

if (lambda ==0) {
  beta.star <- SC.post
} else{
  
  g <- rbind(diag(x=1,N,N),Sigma.true, -Sigma.true)
  
  h <- rbind(matrix(0,N,1),cbind(c(Sigma.true%*%SC.true-lambda,
                                   -(Sigma.true%*%SC.true+lambda))))
  beta.star <- limSolve::lsei(A=mu1,B=tau+mu1%*%SC.post,E=matrix(1,1,N),F=1,
                              G=g,H=h,type=2)$X
  
  # muY <- as.numeric(tau+mu1%*%SC.post)
  # betaHat <- CVXR::Variable(N)
  # 
  # objective   <- CVXR::Minimize( CVXR::quad_form(betaHat,mu1%*%t(mu1)) -2*sum(muY*mu1*betaHat))
  # constraints <- list(betaHat>=0,sum(betaHat)==1, sum(abs(Sigma.true%*%(betaHat-SC.true))^2)<=lambda^2)
  # problem     <- CVXR::Problem(objective, constraints)
  # skip_to_next <- F
  # suppressWarnings({
  #   solution <- CVXR::solve(problem)
  #   tryCatch( beta.cvxr <- as.vector(solution$getValue(betaHat)),
  #             error = function(e) { skip_to_next <<- TRUE})
  #   if(skip_to_next) { next }
  #   # if(any(is.na(beta.m))) {next}
  # })
  
  
}
# tau.star <- tau+mu1%*%(SC.post-beta.cvxr)
tau.star <-c(tau+mu1%*%(SC.post-beta.star))
# tau.star
# if (setting ==1) {
#   lambda <- 0
#   # tau.star <- tau
# }
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
  sig <- sqrt(mean((sc(Y0,X0)$u.hat)^2))
  tau.SC <- mean(Y1-X1%*%SC.beta)
  const.dep <- ifelse(phi==0.5,1,1/2)
  # CI.SC <- c(tau.SC-qnorm(1-alpha/2)*sig/sqrt(T1),tau.SC+qnorm(1-alpha/2)*sig/sqrt(T1))
  ARSC.fit <- ARSC(Y0,Y1,X0,X1,lambda = lambda, method = 'specify',dependent = dependent,c = const.dep,
                   Inference = T,true_beta = beta.star,true_mu = mu1,M = 500)
  ARSC.fit
  
  # CI.SC
  
  # ARSC(Y0,Y1,X0,X1,lambda = lambda, method = 'specify',nu = 0.1*(1+lambda),Inference = T)
  # ARSC.fit$tauHat
  # ARSC.fit2 <- ARSC(Y0,Y1,X0,X1,lambda = lambda, eta = 0.1*sqrt(N*log(N)/T0),method = 'specify',nu = 1e-2)
  dt <- rbind(dt,data.frame(setting = setting,T0 = T0,T1 = T1,tauHat = ARSC.fit$tauHat,
                            tau = tau,tau.star = tau.star,check = ARSC.fit$check, method = 'ARSC'),
              data.frame(setting = setting,T0 = T0,T1 = T1,tauHat = tau.SC,
                         tau = tau,tau.star = tau.star,check = 1, method = 'SC'))
  dt.CI <- rbind(dt.CI,
                 data.frame(setting = setting,T0 = T0,T1 = T1,CI.tau = rbind(ARSC.fit$CI.tau), 
                            negli.m = ARSC.fit$negli.min, 
                            negli.est = abs(sum(-(mu1hat - mu1)*beta.star))+abs(sum(mu1hat*(beta.star - ARSC.fit$betaHat))),
                            tau = tau,tau.star = tau.star,method = 'ARSC'))
  # l2.beta[m,] <- ARSC.fit$beta.min.2
  # linf.beta[m,] <- ARSC.fit$beta.min.inf
  # mu.beta[m,] <- ARSC.fit$beta.min.mu
  # l2.vec[m] <- ARSC.fit$l2.min
  # linf.vec[m] <- ARSC.fit$linf.min
  # mu.vec[m] <- ARSC.fit$mu.min
  # beta.mat[m,] <- ARSC.fit$betaHat
}

filename <- paste("setting",setting,"-phi",phi,"-pre",T0,"-post",T1,"-tau",tau,".RData", sep="")
setwd("~/ARSC/resource")
save.image(filename)

sqrt(max(eigen(t(X0)%*%X0/T0)$values)/max(min(eigen(t(X0)%*%X0/T0)$values),1e-3))
sqrt(max(eigen(Sigma.true)$values)/max(min(eigen(Sigma.true)$values),1e-3))

############################
# Estimator #
tau.star
# setwd("~/Dropbox/Taehyeon/DRO/R-code")
# save.image(file = "sim-dt.Rdata")
# library(ggplot2)
appender_t1 <- function(string) latex2exp::TeX(paste("$T_1 = $", string)) 
appender_t0 <- function(string) latex2exp::TeX(paste("$T_0 = $", string)) 

setting <- 1
t <--.5
t00 <- 200
# mean(dt[dt$tau ==t&dt$T0==t00&dt$setting == setting&dt$method == 'ARSC'&dt$T1==200,]$tauHat)
# unique(dt[dt$tau ==t&dt$T0==t00&dt$setting == setting&dt$method == 'ARSC'&dt$T1==200,]$tau.star)
ggplot(dt[dt$tau ==t&dt$T0==t00&dt$setting == setting,],
       aes(x=method, y=tauHat,fill = method)) +
  # scale_y_discrete()+
  geom_violin(alpha = 0.5,scale = 'width')+
  geom_hline(aes(yintercept = tau.star),colour ='red',linetype="dashed")+
  geom_hline(aes(yintercept = tau),colour ='blue',linetype="dashed")+
  ylab(latex2exp::TeX("$\\hat{\\tau}$"))+
  facet_wrap(~T1,
             labeller = as_labeller(appender_t1,
                                    default = label_parsed),
             ncol = 4)+
  theme(text = element_text(size = 20)) + guides(fill="none")+
  ggtitle(latex2exp::TeX(paste0("(S",setting,"), T0 = ",t00,", \\tau = ",t)))+
  theme(plot.title = element_text(hjust = 0.5))


# How many times increase threshold
ggplot(dt[dt$tau ==t&dt$T0==t00&dt$setting == setting&dt$method == 'ARSC',],
       aes(x=check)) +
  geom_histogram()+
  # scale_y_discrete()+
  # geom_violin(alpha = 0.5,scale = 'width',trim = F)+
  facet_wrap(~T1,
             labeller = as_labeller(appender_t1,
                                    default = label_parsed),
             ncol = 4)+
  theme(text = element_text(size = 20)) + guides(fill="none")+
  ggtitle(latex2exp::TeX(paste0("(S",setting,"), T0 = ",t00,", \\tau = ",t)))+
  theme(plot.title = element_text(hjust = 0.5))



ggplot(dt[dt$tau ==t&dt$T1==t00&dt$setting == setting,],
       aes(x=method, y=tauHat,fill = method)) +
  # scale_y_discrete()+
  geom_violin(alpha = 0.5,scale = 'width',trim = F)+
  geom_hline(aes(yintercept = tau.star),colour ='red',linetype="dashed")+
  geom_hline(aes(yintercept = tau),colour ='blue',linetype="dashed")+
  ylab(latex2exp::TeX("$\\hat{\\tau}$"))+
  facet_wrap(~T0,
             labeller = as_labeller(appender_t0,
                                    default = label_parsed),
             ncol = 4)+
  theme(text = element_text(size = 20)) + guides(fill="none")+
  ggtitle(latex2exp::TeX(paste0("(S",setting,"), T1 = ",t00,", \\tau = ",t)))+
  theme(plot.title = element_text(hjust = 0.5))
setwd("/Users/taehyeon/Library/CloudStorage/Dropbox/Taehyeon/DRO/Figures/Estimators")


for (setting in 1:3) {
  for (t00 in c(50,100,200)) {
    for (t in c(-1,-0.5,0.5,1)) {
      temp <- ggplot(dt[dt$setting ==setting
                              & dt$T0 ==t00
                              & dt$method %in% c("ARSC","SC")
                              & dt$T1 %in% c(50,100,200)
                              & dt$tau ==t,],
                     aes(x=method, y=tauHat,fill = method)) +
        # scale_y_discrete()+
        geom_violin(alpha = 0.5,scale = 'width',trim = F)+
        geom_hline(aes(yintercept = tau.star),colour ='red',linetype="dashed")+
        geom_hline(aes(yintercept = tau),colour ='blue',linetype="dashed")+
        ylab(latex2exp::TeX("$\\hat{\\tau}$"))+
        facet_wrap(~T1,
                   labeller = as_labeller(appender_t1, 
                                          default = label_parsed),
                   ncol = 3)+
        theme(text = element_text(size = 20)) + guides(fill="none")+
        xlab("")+
        ggtitle(latex2exp::TeX(paste0("(S",setting,"), T0 = ",t00,", \\tau = ",t)))+
        theme(plot.title = element_text(hjust = 0.5))
      # print(temp)
      plotname <- paste0("setting",setting,"-tau",t,"-pre",t00,".png")
      ggsave(plot = temp,filename = plotname)
    }
  }
}


####################
# Inference #
################
t0 <- 200
t1 <- 200
s <- 3
t <- 1
SC.idx.oracle <- dt$setting==s&dt$T0==t0&dt$T1==t1&dt$tau==t&dt$method=="SC"
# SC.idx <- dt.CI$setting==s&dt.CI$T0==t0&dt.CI$T1==t1&dt.CI$tau==t&dt.CI$method=="SC"
ARSC.idx.oracle <-  dt$setting==s&dt$T0==t0&dt$T1==t1&dt$tau==t&dt$method=="ARSC"
ARSC.idx <- dt.CI$setting==s&dt.CI$T0==t0&dt.CI$T1==t1&dt.CI$tau==t&dt.CI$method=="ARSC" # truncated
ARSC.idx2 <- dt.CI$setting==s&dt.CI$T0==t0&dt.CI$T1==t1&dt.CI$tau==t&dt.CI$method=="ARSC2" 
chi <- qchisq(1-alpha, 1, ncp = unique(((dt[SC.idx.oracle,]$tauHat-dt[SC.idx.oracle,]$tau)/sd(dt[SC.idx.oracle,]$tauHat))^2), 
              lower.tail = TRUE, log.p = FALSE)
chi.ARSC <- qchisq(1-alpha, 1, ncp = unique(((dt[ARSC.idx.oracle,]$tauHat-dt[ARSC.idx.oracle,]$tau.star)/sd(dt[ARSC.idx.oracle,]$tauHat))^2), 
                       lower.tail = TRUE, log.p = FALSE)
SC.OBA<- cbind(dt[SC.idx.oracle,]$tauHat-sd(dt[SC.idx.oracle,]$tauHat)*sqrt(chi),
               dt[SC.idx.oracle,]$tauHat+sd(dt[SC.idx.oracle,]$tauHat)*sqrt(chi)) # OBA
ARSC.OBA <- cbind(dt[ARSC.idx.oracle,]$tauHat-sd(dt[ARSC.idx.oracle,]$tauHat)*sqrt(chi.ARSC),
                  dt[ARSC.idx.oracle,]$tauHat+sd(dt[ARSC.idx.oracle,]$tauHat)*sqrt(chi.ARSC)) # OBA
# Coverage #

# SC estimator

# ARSC

# ARSC 2
# Oracle


# Length ratio



# Length ratio





for (s in 1:3) {
  for (t in c(-1,-0.5,0,0.5,1)) {
    cat("
\\begin{table}[]
\\centering
\\resizebox{.9\\textwidth}{!}{%
\\begin{tabular}{cc|cccc|cccc}
\\multicolumn{2}{c|}{\\multirow{2}{*}{",paste0("(S",s,"), $\\tau = ",t,"$"),"}} & \\multicolumn{4}{c|}{Coverage} & \\multicolumn{4}{c}{Length} \\\\ \\cline{3-10}
\\multicolumn{2}{c|}{}                                        & OBA($\\bar{\\tau}$)     & OBA($\\tau^*$)   & CI  & Int &  OBA($\\bar{\\tau}$)  & OBA($\\tau^*$)&  CI& Int \\\\ \\hline")
    for (t0 in c(50,100,200)) {
      for (t1 in c(50,100,200)) {
        SC.idx.oracle <- dt$setting==s&dt$T0==t0&dt$T1==t1&dt$tau==t&dt$method=="SC"
        # SC.idx <- dt.CI$setting==s&dt.CI$T0==t0&dt.CI$T1==t1&dt.CI$tau==t&dt.CI$method=="SC"
        ARSC.idx.oracle <-  dt$setting==s&dt$T0==t0&dt$T1==t1&dt$tau==t&dt$method=="ARSC"
        
        ARSC.idx <- dt.CI$setting==s&dt.CI$T0==t0&dt.CI$T1==t1&dt.CI$tau==t&dt.CI$method=="ARSC" # truncated
        ARSC.idx2 <- dt.CI$setting==s&dt.CI$T0==t0&dt.CI$T1==t1&dt.CI$tau==t&dt.CI$method=="ARSC2" 
        chi <- qchisq(1-alpha, 1, 
                      ncp = unique(((mean(dt[SC.idx.oracle,]$tauHat)-dt[SC.idx.oracle,]$tau)/sd(dt[SC.idx.oracle,]$tauHat))^2), 
                      lower.tail = TRUE, log.p = FALSE)
        chi.ARSC <- qchisq(1-alpha, 1, ncp = unique(((mean(dt[ARSC.idx.oracle,]$tauHat)-dt[ARSC.idx.oracle,]$tau.star)/sd(dt[ARSC.idx.oracle,]$tauHat))^2), 
                           lower.tail = TRUE, log.p = FALSE)
        SC.OBA <- cbind(dt[SC.idx.oracle,]$tauHat-sd(dt[SC.idx.oracle,]$tauHat)*sqrt(chi),
                           dt[SC.idx.oracle,]$tauHat+sd(dt[SC.idx.oracle,]$tauHat)*sqrt(chi)) # OBA
        ARSC.OBA <- cbind(dt[ARSC.idx.oracle,]$tauHat-sd(dt[ARSC.idx.oracle,]$tauHat)*sqrt(chi.ARSC),
                          dt[ARSC.idx.oracle,]$tauHat+sd(dt[ARSC.idx.oracle,]$tauHat)*sqrt(chi.ARSC)) # OBA
        tt <- unique(dt[ARSC.idx.oracle,]$tau.star)
        if (t1==50 & t0==50) {
          cat("\\multicolumn{1}{c|}{\\multirow{3}{*}{$T_0 = 50$}}  & $T_1 = 50$  & ",
              sum(apply(SC.OBA,1,function(x){ifelse(x[1]<t & x[2]>t,1,0)}))/nsims,"&",
              sum(apply(ARSC.OBA,1,function(x){ifelse(x[1]<tt & x[2]>tt,1,0)}))/nsims,"&",
              length(which(dt.CI[ARSC.idx,]$CI.tau.1<=round(unique(dt.CI[ARSC.idx,]$tau.star),3) &
                             dt.CI[ARSC.idx,]$CI.tau.2>=round(round(unique(dt.CI[ARSC.idx,]$tau.star),3),3)))/nsims,"&",
              length(which(dt.CI[ARSC.idx2,]$CI.tau.1<=unique(dt.CI[ARSC.idx2,]$tau.star) &
                             dt.CI[ARSC.idx2,]$CI.tau.2>=round(unique(dt.CI[ARSC.idx2,]$tau.star),3)))/nsims,"&",
              round(mean(apply(SC.OBA,1,diff)),3),"&",
              round(mean(apply(ARSC.OBA,1,diff)),3),"&",
              round(mean(dt.CI[ARSC.idx,]$CI.tau.2-dt.CI[ARSC.idx,]$CI.tau.1),3),"&",
              
              round(mean(dt.CI[ARSC.idx2,]$CI.tau.2-dt.CI[ARSC.idx2,]$CI.tau.1),3)
              ," \\\\ ")
        }
        else if (t0==50 & t1==100) {
          cat("\\multicolumn{1}{c|}{}                             & $T_1 = 100$ & ",
              sum(apply(SC.OBA,1,function(x){ifelse(x[1]<t & x[2]>t,1,0)}))/nsims,"&",
              
              sum(apply(ARSC.OBA,1,function(x){ifelse(x[1]<tt & x[2]>tt,1,0)}))/nsims,"&",
              length(which(dt.CI[ARSC.idx,]$CI.tau.1<=round(unique(dt.CI[ARSC.idx,]$tau.star),3) &
                             dt.CI[ARSC.idx,]$CI.tau.2>=round(round(unique(dt.CI[ARSC.idx,]$tau.star),3),3)))/nsims,"&",
              length(which(dt.CI[ARSC.idx2,]$CI.tau.1<=unique(dt.CI[ARSC.idx2,]$tau.star) &
                             dt.CI[ARSC.idx2,]$CI.tau.2>=round(unique(dt.CI[ARSC.idx2,]$tau.star),3)))/nsims,"&",
              round(mean(apply(SC.OBA,1,diff)),3),"&",
              round(mean(apply(ARSC.OBA,1,diff)),3),"&",
              round(mean(dt.CI[ARSC.idx,]$CI.tau.2-dt.CI[ARSC.idx,]$CI.tau.1),3),"&",
              
              round(mean(dt.CI[ARSC.idx2,]$CI.tau.2-dt.CI[ARSC.idx2,]$CI.tau.1),3)," \\\\ ")
        }
        else if (t1==200 & t0!=200){
          cat("\\multicolumn{1}{c|}{}                             & $T_1 = 200$ &",
              sum(apply(SC.OBA,1,function(x){ifelse(x[1]<t & x[2]>t,1,0)}))/nsims,"&",
              
              sum(apply(ARSC.OBA,1,function(x){ifelse(x[1]<tt & x[2]>tt,1,0)}))/nsims,"&",
              length(which(dt.CI[ARSC.idx,]$CI.tau.1<=round(unique(dt.CI[ARSC.idx,]$tau.star),3) &
                             dt.CI[ARSC.idx,]$CI.tau.2>=round(round(unique(dt.CI[ARSC.idx,]$tau.star),3),3)))/nsims,"&",
              length(which(dt.CI[ARSC.idx2,]$CI.tau.1<=unique(dt.CI[ARSC.idx2,]$tau.star) &
                             dt.CI[ARSC.idx2,]$CI.tau.2>=round(unique(dt.CI[ARSC.idx2,]$tau.star),3)))/nsims,"&",
              round(mean(apply(SC.OBA,1,diff)),3),"&",
              round(mean(apply(ARSC.OBA,1,diff)),3),"&",
              round(mean(dt.CI[ARSC.idx,]$CI.tau.2-dt.CI[ARSC.idx,]$CI.tau.1),3),"&",
              
              round(mean(dt.CI[ARSC.idx2,]$CI.tau.2-dt.CI[ARSC.idx2,]$CI.tau.1),3)," \\\\ \\hline")
        }
        else if (t1==50 & t0==100) {
          cat("\\multicolumn{1}{c|}{\\multirow{3}{*}{$T_0 = 100$}} & $T_1 = 50$  &",
              sum(apply(SC.OBA,1,function(x){ifelse(x[1]<t & x[2]>t,1,0)}))/nsims,"&",
              
              sum(apply(ARSC.OBA,1,function(x){ifelse(x[1]<tt & x[2]>tt,1,0)}))/nsims,"&",
              length(which(dt.CI[ARSC.idx,]$CI.tau.1<=round(unique(dt.CI[ARSC.idx,]$tau.star),3) &
                             dt.CI[ARSC.idx,]$CI.tau.2>=round(round(unique(dt.CI[ARSC.idx,]$tau.star),3),3)))/nsims,"&",
              length(which(dt.CI[ARSC.idx2,]$CI.tau.1<=unique(dt.CI[ARSC.idx2,]$tau.star) &
                             dt.CI[ARSC.idx2,]$CI.tau.2>=round(unique(dt.CI[ARSC.idx2,]$tau.star),3)))/nsims,"&",
              round(mean(apply(SC.OBA,1,diff)),3),"&",
              round(mean(apply(ARSC.OBA,1,diff)),3),"&",
              round(mean(dt.CI[ARSC.idx,]$CI.tau.2-dt.CI[ARSC.idx,]$CI.tau.1),3),"&",
              
              round(mean(dt.CI[ARSC.idx2,]$CI.tau.2-dt.CI[ARSC.idx2,]$CI.tau.1),3)," \\\\
")
        }
        else if (t1==50 & t0==200) {
          cat("\\multicolumn{1}{c|}{\\multirow{3}{*}{$T_0 = 200$}} & $T_1 = 50$  &", 
              sum(apply(SC.OBA,1,function(x){ifelse(x[1]<t & x[2]>t,1,0)}))/nsims,"&",
              
              sum(apply(ARSC.OBA,1,function(x){ifelse(x[1]<tt & x[2]>tt,1,0)}))/nsims,"&",
              length(which(dt.CI[ARSC.idx,]$CI.tau.1<=round(unique(dt.CI[ARSC.idx,]$tau.star),3) &
                             dt.CI[ARSC.idx,]$CI.tau.2>=round(round(unique(dt.CI[ARSC.idx,]$tau.star),3),3)))/nsims,"&",
              length(which(dt.CI[ARSC.idx2,]$CI.tau.1<=unique(dt.CI[ARSC.idx2,]$tau.star) &
                             dt.CI[ARSC.idx2,]$CI.tau.2>=round(unique(dt.CI[ARSC.idx2,]$tau.star),3)))/nsims,"&",
              round(mean(apply(SC.OBA,1,diff)),3),"&",
              round(mean(apply(ARSC.OBA,1,diff)),3),"&",
              round(mean(dt.CI[ARSC.idx,]$CI.tau.2-dt.CI[ARSC.idx,]$CI.tau.1),3),"&",
              
              round(mean(dt.CI[ARSC.idx2,]$CI.tau.2-dt.CI[ARSC.idx2,]$CI.tau.1),3)," \\\\
")
        }
        else if ( t0==200 & t1==200){
          cat("\\multicolumn{1}{c|}{}                             & $T_1 = 200$ &",
              sum(apply(SC.OBA,1,function(x){ifelse(x[1]<t & x[2]>t,1,0)}))/nsims,"&",
              
              sum(apply(ARSC.OBA,1,function(x){ifelse(x[1]<tt & x[2]>tt,1,0)}))/nsims,"&",
              length(which(dt.CI[ARSC.idx,]$CI.tau.1<=round(unique(dt.CI[ARSC.idx,]$tau.star),3) &
                             dt.CI[ARSC.idx,]$CI.tau.2>=round(round(unique(dt.CI[ARSC.idx,]$tau.star),3),3)))/nsims,"&",
              length(which(dt.CI[ARSC.idx2,]$CI.tau.1<=unique(dt.CI[ARSC.idx2,]$tau.star) &
                             dt.CI[ARSC.idx2,]$CI.tau.2>=round(unique(dt.CI[ARSC.idx2,]$tau.star),3)))/nsims,"&",
              round(mean(apply(SC.OBA,1,diff)),3),"&",
              round(mean(apply(ARSC.OBA,1,diff)),3),"&",
              round(mean(dt.CI[ARSC.idx,]$CI.tau.2-dt.CI[ARSC.idx,]$CI.tau.1),3),"&",
              
              round(mean(dt.CI[ARSC.idx2,]$CI.tau.2-dt.CI[ARSC.idx2,]$CI.tau.1),3),"          
")
        }
        else{
          cat("\\multicolumn{1}{c|}{}                             & $T_1 = 100$ &",
              sum(apply(SC.OBA,1,function(x){ifelse(x[1]<t & x[2]>t,1,0)}))/nsims,"&",
              
              sum(apply(ARSC.OBA,1,function(x){ifelse(x[1]<tt & x[2]>tt,1,0)}))/nsims,"&",
              length(which(dt.CI[ARSC.idx,]$CI.tau.1<=round(unique(dt.CI[ARSC.idx,]$tau.star),3) &
                             dt.CI[ARSC.idx,]$CI.tau.2>=round(round(unique(dt.CI[ARSC.idx,]$tau.star),3),3)))/nsims,"&",
              length(which(dt.CI[ARSC.idx2,]$CI.tau.1<=unique(dt.CI[ARSC.idx2,]$tau.star) &
                             dt.CI[ARSC.idx2,]$CI.tau.2>=round(unique(dt.CI[ARSC.idx2,]$tau.star),3)))/nsims,"&",
              round(mean(apply(SC.OBA,1,diff)),3),"&",
              round(mean(apply(ARSC.OBA,1,diff)),3),"&",
              round(mean(dt.CI[ARSC.idx,]$CI.tau.2-dt.CI[ARSC.idx,]$CI.tau.1),3),"&",
              
              round(mean(dt.CI[ARSC.idx2,]$CI.tau.2-dt.CI[ARSC.idx2,]$CI.tau.1),3)," \\\\
")
        }
        
        
      }
    }
    cat("\\end{tabular}%
              }
              \\caption{",paste0("Inference result for (S",s,"), $\\tau = ",t,"$"),"}
              \\label{tab:",paste0("S",s," tau",t),"}
            \\end{table}\n\n")

  
  
  }
}


setwd("~/Dropbox/Taehyeon/DRO/R-code")
save.image(file = "sim-dt.Rdata")

temp <- generate_AR1_process(T0+T1, 1, .5, 1, 0)
var(temp)
sandwich::NeweyWest(lm(temp~1),lag = 10)
sqrt(as.numeric(HAC_cov_matrix(temp)))

