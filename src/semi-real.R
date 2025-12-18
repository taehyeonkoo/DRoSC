
######################
# Basque

library(Synth)
library(ggcorrplot)
source("~/Dropbox/Taehyeon/Synthetic Control/R-code/helpers.R")
### Basque ###
### T0 = 15, N = 16, T1 = 28
data("basque") 

T0 <- 15
T1 <- 28
Y <- basque[basque$regionname=="Basque Country (Pais Vasco)",]$gdpcap
X <- matrix(basque[basque$regionname!="Basque Country (Pais Vasco)" & basque$regionname!="Spain (Espana)",]$gdpcap,
            nrow = T0+T1,byrow = F)
colnames(X) <- unique(basque[basque$regionname!="Basque Country (Pais Vasco)"& basque$regionname!="Spain (Espana)",]$regionname)
colnames(X) <-c("Andalucia","Aragon","Asturias","Baleares","Canarias","Cantabria",
                "Castilla Y Leon","Castilla-La Mancha","Cataluna","Comunidad Valenciana","Extremadura","Galicia",
                "Madrid","Murcia","Navarra","Rioja"
                )
Y0 <- Y[1:T0]
Y1 <- Y[-(1:T0)]
X0 <- X[1:T0,]
X1 <- X[-(1:T0),]

N <- ncol(X)
SC.fit <- sc(Y0,X0)
SC.beta <- SC.fit$w.hat

# Correlation plot (Figure 1)#
SC.names <- names(SC.beta[which(SC.beta>0)][order(SC.beta[which(SC.beta>0)],decreasing = F)])
cor.plot <- ggcorrplot(cor(X0)[,SC.names])+ theme( legend.text = element_text(size = 14),text = element_text(size = 15))
cor.plot

# Adding noise the pre-treatment data #
c.cand <- c(0,0.05,0.1,0.15)
dt.var <- dt <- data.frame()
nsims <- 1000
beta.mat.DRSC <- beta.mat <- matrix(nrow = nsims,ncol = N)
beta.list.DRSC <- beta.list <- list()
for (c in c.cand) {
  cat("c = ",c,"\n")
  # var.mat <- matrix(nrow = length(l.cand),ncol = 2)
  sc.vec <- vector()
  drsc.vec <- vector()
  for (nsim in 1:nsims) {
    if (nsim %% (nsims/2) == 1) {
      cat("nsim = ",nsim,"\n")
    }
    X0.new <- X0+MASS::mvrnorm(T0,mu = rep(0,N), Sigma = c*diag(sqrt(diag(cov(X0)))))
    Y0.new <- Y0+rnorm(T0,sd = c*sd(Y0))
    beta.mat[nsim,] <- sc(Y0.new,X0.new)$w.hat
    tau.sc <- mean(Y1-X1%*%sc(Y0.new,X0.new)$w.hat)
    sc.vec[nsim] <- tau.sc
    ARSC.fit <- ARSC(Y0.new,Y1,X0.new,X1,lambda = 0)
    drsc.vec[nsim] <- ARSC.fit$tauHat
    beta.mat.DRSC[nsim,] <- ARSC.fit$betaHat
  }
  beta.list[[which(c.cand==c)]] <- beta.mat
  beta.list.DRSC[[which(c.cand==c)]] <- beta.mat.DRSC
  dt <- rbind(dt,data.frame(tauHat = drsc.vec, c = c, method = 'DRoSC'),
              data.frame(tauHat = sc.vec, c = c,method = 'SC'))
}


# Columns to consider based on SC.beta
cols_to_keep <- which(apply(beta.list.DRSC[[2]],2,function(col) mean(col > 0)>0))
cols_to_keep <-cols_to_keep[order(SC.beta[cols_to_keep], decreasing = TRUE)]
cols_to_keep <- setdiff(cols_to_keep,which(names(SC.beta)%in%c("Navarra","Canarias","Murcia")))
rows <- list()

# Loop over each threshold
for (j in seq_along(c.cand)) {
  c_val <- c.cand[j]
  
  # Compute the proportion for each region (column)
  prop <- apply(beta.list.DRSC[[j]][, cols_to_keep, drop = FALSE], 2, function(x) mean(x > c_val))
  
  # Store the result with row name as threshold
  rows[[j]] <- prop
}

# Combine the list of proportions into a data frame
result_df <- do.call(rbind, rows)
result_df <- as.data.frame(result_df)
# Assign column names (regions)
colnames(result_df) <- names(SC.beta)[cols_to_keep]

# Add the threshold column by directly assigning it as a new column
result_df$threshold <- factor(c.cand)

# If you want to use the long format for plotting
df_long <- result_df %>%
  pivot_longer(-threshold, names_to = "region", values_to = "proportion")
df_long$region <- factor(df_long$region, levels = names(SC.beta)[cols_to_keep])
# Plot: x = threshold, y = proportion, facet by region
plot.basque.supp <- ggplot(df_long, aes(x = threshold, y = proportion, fill = region)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ region, nrow = 1) +
  labs(x = "c", y = "Proportion") +
  theme_minimal()+
  theme(legend.position = "none",text = element_text(size = 15))+
  theme(strip.text = element_text(size = 20))
# Figure 2 #
plot.basque.supp


### For DRoSC ###
# Columns to consider based on SC.beta
cols_to_keep <- which(apply(beta.list.DRSC[[1]],2,function(col) mean(col > 0)>0.05))#which(SC.beta>0)
cols_to_keep <-cols_to_keep[order(SC.beta[cols_to_keep], decreasing = TRUE)]
cols_to_keep <- setdiff(cols_to_keep,which(names(SC.beta)%in%c("Navarra","Canarias","Murcia")))
rows <- list()

# Loop over each threshold
for (j in seq_along(c.cand)) {
  c_val <- c.cand[j]
  
  # Compute the proportion for each region (column)
  prop <- apply(beta.list.DRSC[[j]][, cols_to_keep, drop = FALSE], 2, function(x) mean(x > c_val))
  
  # Store the result with row name as threshold
  rows[[j]] <- prop
}

# Combine the list of proportions into a data frame
result_df <- do.call(rbind, rows)
X0[,cols_to_keep]
result_df <- as.data.frame(result_df)
# Assign column names (regions)
colnames(result_df) <- names(SC.beta)[cols_to_keep]

# Add the threshold column by directly assigning it as a new column
result_df$threshold <- factor(c.cand)

# If you want to use the long format for plotting
df_long <- result_df %>%
  pivot_longer(-threshold, names_to = "region", values_to = "proportion")
df_long$region <- factor(df_long$region, levels = names(SC.beta)[cols_to_keep])
# Plot: x = threshold, y = proportion, facet by region
plot.basque.supp <- ggplot(df_long, aes(x = threshold, y = proportion, fill = region)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ region, nrow = 1) +
  labs(x = "c", y = "Proportion") +
  theme_minimal()+
  theme(legend.position = "none",text = element_text(size = 15))+
  theme(strip.text = element_text(size = 20))
ggsave(plot = plot.basque.supp, 
       filename = "Basque_supp.png",
       device = "png", width = 5, height = 2)


# Weight shift #
# g <- 0.1
set.seed(42)
mu0 <- apply(X0,2,mean)
mu1 <- apply(X1,2,mean)
prop.dt <- data.frame()
nsims <- 1000
k.cand <- c(0.05,0.1,0.2,0.3,0.4)#seq(0.05,0.25,0.05)
for (k in k.cand) {
# k <- 0.05
  SC.post <- SC.beta
  SC.post[which(colnames(X0)%in%c("Cataluna"))] <- k*SC.beta[which(colnames(X0)%in%c("Baleares"))]
  SC.post[which(colnames(X0)%in%c("Asturias"))] <-k*SC.beta[which(colnames(X0)%in%c("Rioja"))]
  
  SC.post[which(colnames(X0)%in%c("Baleares"))] <- (1-k)*SC.beta[which(colnames(X0)%in%c("Baleares"))]
  SC.post[which(colnames(X0)%in%c("Rioja"))] <- (1-k)*SC.beta[which(colnames(X0)%in%c("Rioja"))]
  
  stopifnot(sum(SC.post)|all(SC.post>=0) )
   # sum-to-one
  # all non-negative
  tau.t <- Y1-X1%*%SC.post
  tau.SC<- mean(Y1-X1%*%SC.post)
  pre.var <- mean((Y0-X0%*%SC.beta)^2)
  # post.var <- mean((Y1-tau.SC-X1%*%SC.beta)^2)
  
  
  X0.var <- min(diag(cov(X0)))
  X1.var <- min(diag(cov(X1)))
  
  arsc.vec <- tau.vec <- vector()
  Sigma.true <- t(X0)%*%X0/T0+diag(X0.var,N) # Covariance + mu
  # eigen(Sigma.true)$values
  lambda <- max(abs((Sigma.true%*%(SC.beta-SC.post))))
  
  lambda
  
  g <- rbind(diag(x=1,N,N),Sigma.true, -Sigma.true)
  h <- rbind(matrix(0,N,1),cbind(c(Sigma.true%*%SC.beta-lambda,
                                   -(Sigma.true%*%SC.beta+lambda))))
  
  if (all(SC.beta==SC.post)) {
    beta.star <- SC.post
  } else{
    beta.star <- limSolve::lsei(A=mu1,B=tau.SC+mu1%*%SC.post,E=matrix(1,1,N),F=1,
                                G=g,H=h,type=2)$X
  }
  
  
  # SC.post
  # beta.star
  tau.star <-c(tau.SC+mu1%*%(SC.post-beta.star)) # tau star in semi-real data analysis
  # tau.star
  # print(tau.star)
  max(abs(Sigma.true%*%(beta.star-SC.beta)))
  
  cov(X0)
  
  # semi.tau <- MASS::mvrnorm(nsims,mu = tau.t, Sigma = diag(as.numeric(var(tau.t)),T1))
  Y0.mat <- matrix(nrow = nsims, ncol = T0)
  Y1.mat <- matrix(nrow = nsims, ncol = T1)
  Y1.obs.mat <- matrix(nrow = nsims, ncol = T1)
  SC.mat <- matrix(nrow = nsims,ncol = N)
  CI.SC.mat <- matrix(nrow = nsims,ncol = 2)
  CI.ARSC.mat<- matrix(nrow = nsims,ncol = 2)
  length.SC <- rep(0,nsims)
  length.ARSC <- rep(0,nsims)
  cnt.SC <- 0
  cnt.ARSC <- 0
  CI.mat <- matrix(nrow = nsims, ncol = 2)
  beta.mat <- matrix(nrow = nsims,ncol = length(beta.star))
  M <- 500
  
  
  for (i in 1:nsims) {
    semi.X0 <- X0+matrix(rnorm(T0*N,mean = 0, sd = sqrt(X0.var)),
                         nrow = T0,ncol = N)
    semi.Y0 <- as.vector(semi.X0%*%SC.beta+rnorm(T0,sd = sqrt(pre.var)))
    # Y0.mat[i,] <- semi.Y0
    # for (t in 1:T1) {
    #   semi.X1[t,] <- MASS::mvrnorm(n = 1,X1[t,],new.X1)
    # }
    semi.X1 <- X1+matrix(rnorm(T1*N,mean = 0, sd = sqrt(X1.var)),
                         nrow = T1,ncol = N)
    semi.Y1 <- semi.X1%*%SC.post+rnorm(T1,sd = sqrt(pre.var))
    # Y1.mat[i,] <- semi.Y1
    delta.t <- MASS::mvrnorm(1,mu = rep(0,T1), Sigma = diag(as.numeric(var(tau.t)),T1))
    semi.Y1.obs <- as.vector(semi.Y1+tau.t+delta.t)
    # Y1.obs.mat[i,] <- semi.Y1.obs
    SC.i <- sc(semi.Y0,semi.X0)$w.hat
    sig <- sd(sc(semi.Y0,semi.X0)$u.hat)
    # SC.i
    # semi.tau.t <- rnorm(T1,tau.SC,sd = sd(tau.t))
    # SC.mat[i,] <- SC.i
    tau.i <- mean(semi.Y1.obs-semi.X1%*%SC.i)
    tau.vec[i] <- tau.i
    arsc.i <- ARSC(semi.Y0,semi.Y1.obs,semi.X0,semi.X1,lambda = lambda,
                   method = 'specify',Inference = F,M = M,
                   true_beta = beta.star, true_mu = apply(X1,2,mean))

    arsc.vec[i] <- arsc.i$tauHat
    
  }
  prop.dt <- rbind(prop.dt, data.frame(tauHat = tau.vec, estimand = tau.SC, k = k, method = 'SC'),
        data.frame(tauHat = arsc.vec, estimand = tau.star, k = k, method = 'DRoSC'))
 
}
library(dplyr)
# Filter for SC method only
sc_only <- prop.dt %>% filter(method == "SC")

# Plot violin plot of tauHat by k
ggplot(sc_only, aes(x = factor(k), y = tauHat)) +
  geom_violin(fill = "skyblue") +
  labs(x = expression(kappa), y = expression(hat(tau)^{SC})) +
  # geom_hline(aes(yintercept = estimand), color = "blue", linetype = "dashed") +
  geom_segment(aes(x = as.numeric(factor(k)) - 0.5,
                   xend = as.numeric(factor(k)) + 0.5,
                   y = estimand,
                   yend = estimand),
               color = "blue",  linetype = "dashed")+
  theme_minimal()+
    theme(text = element_text(size = 15))

# Filter for SC method only
arsc_only <- prop.dt %>% filter(method == "DRoSC")

# Plot violin plot of tauHat by k
ggplot(arsc_only, aes(x = factor(k), y = tauHat)) +
  geom_violin(fill = "lightpink",scale = "width") +
  labs(x = expression(kappa), y = expression(hat(tau))) +
  geom_segment(aes(x = as.numeric(factor(k)) - 0.5, 
                   xend = as.numeric(factor(k)) + 0.5,
                   y = estimand, 
                   yend = estimand),
               color = "red",  linetype = "dashed")+
  theme_minimal()+
  theme(text = element_text(size = 15))

  
  


