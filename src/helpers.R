# Generate AR series
# Function to generate N-dimensional AR(1) process
generate_AR1_process <- function(T0, N, phi, Sigma, mu) {
  if (!all(eigen(Sigma)$values > 0)) {
    stop("Sigma must be positive definite.")
  }
  
  # Ensure mu is a vector of length N
  mu <- rep(mu, length.out = N)
  
  # Handle univariate (N=1) and multivariate (N>1) Sigma cases correctly
  if (N == 1) {
    Sigma <- matrix(Sigma, 1, 1)  # Convert scalar to 1x1 matrix
  }
  
  # Compute Cholesky decomposition of Sigma
  # chol_Sigma <- chol(Sigma)
  
  # Initialize the AR(1) process storage
  X <- matrix(0, nrow = T0, ncol = N)
  
  # Generate the initial values for the process
  X[1, ] <- mu + MASS::mvrnorm(1, mu = rep(0, N), Sigma = Sigma)
  
  # Generate AR(1) process for subsequent time points
  for (t in 2:T0) {
    noise <- MASS::mvrnorm(1, mu = rep(0, N), Sigma = Sigma)  # Noise at time t
    X[t, ] <- mu + phi * (X[t - 1, ] - mu) + noise
  }
  return(X)
}
HAC_cov_matrix <- function(Y, lag = NULL, prewhite = F) {
  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)  # Convert to N x 1 matrix
  }
  T0 <- nrow(Y)  # Number of time points
  N <- ncol(Y)  # Dimension of vector process
  Y_bar <- colMeans(Y)  # Sample mean of Y
  
  # Automatic lag selection using Newey-West rule if lag is not provided
  if (is.null(lag)) {
    lag <- floor(4 * (T0 / 100)^(2 / 9))
  }
  
  # Apply prewhitening if needed (AR(1) correction)
  if (prewhite) {
    Y <- Y[-1, , drop = FALSE] - Y[-T0, , drop = FALSE]  # First differences as a simple prewhitening step
    T0 <- nrow(Y)  # Adjust T0 after transformation
  }
  
  # Initialize HAC estimator
  Sigma_HAC <- matrix(0, N, N)
  
  for (h in 0:lag) {
    # Compute autocovariance for lag h
    if (h == 0) {
      Gamma_h <- crossprod(sweep(Y, 2, Y_bar)) / (T0 - 1)  # Variance estimate
    } else {
      Gamma_h <- crossprod(sweep(Y[1:(T0 - h), , drop = FALSE], 2, Y_bar),
                           sweep(Y[(1 + h):T0, , drop = FALSE], 2, Y_bar)) / (T0 - h)
    }
    
    # Bartlett kernel weight
    weight <- 1 - h / (lag + 1)
    
    # Apply weighting
    if (h == 0) {
      Sigma_HAC <- Sigma_HAC + Gamma_h  # No weight for h = 0
    } else {
      Sigma_HAC <- Sigma_HAC + weight * (Gamma_h + t(Gamma_h))  # Symmetric adjustment
    }
  }
  
  return(Sigma_HAC)
}


A1gen <- function(rho,p){
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]<-rho^(abs(i-j))
    }
  }
  return(A1)
}

vecl <- function(A){
  A[upper.tri(A)] <- NA
  A.vec <- as.vector(A)[!is.na(A)]
  return(A.vec)
}


# Standard synthetic control method #
sc <- function(Y0,X0, intercept = FALSE, u = NULL){
  if (length(u)==0) {
    if (intercept) {
      J <- dim(X0)[2]
      X0 <- cbind(1,X0)
      e <- cbind(0,matrix(1,1,J))
      f <- 1
      g <- diag(x=1,J+1,J+1)
      g <- g[-1,]
      h <- matrix(0,J,1)
      w.hat <- limSolve::lsei(A=X0,B=Y0,E=e,F=f,G=g,H=h,type=2)$X
      u.hat <- Y0-X0%*%w.hat
      return(list(u.hat=u.hat,w.hat=w.hat))
    } else {
      J <- dim(X0)[2]
      e <- matrix(1,1,J)
      f <- 1
      g <- diag(x=1,J,J)
      h <- matrix(0,J,1)
      w.hat <- limSolve::lsei(A=X0,B=Y0,E=e,F=f,G=g,H=h,type=2)$X
      u.hat <- Y0-X0%*%w.hat
      return(list(u.hat=u.hat,w.hat=w.hat))
    }
  } else {
    if (intercept) {
      J <- dim(X0)[2]
      X0 <- cbind(1,X0)
      e <- cbind(0,matrix(1,1,J),0)
      f <- 1
      g <- cbind(diag(x=1,J+1,J+1),rep(0,J+1))
      g <- g[-1,]
      h <- matrix(0,J,1)
      w.hat <- limSolve::lsei(A=cbind(X0,u),B=Y0,E=e,F=f,G=g,H=h,type=2)$X
      u.hat <- Y0-X0%*%w.hat
      return(list(u.hat=u.hat,w.hat=w.hat))
    } else{
      J <- dim(X0)[2]
      e <- cbind(matrix(1,1,J),0)
      f <- 1
      g <- cbind(diag(x=1,J,J),rep(0,J))
      h <- rbind(matrix(0,J,1))
      temp <- limSolve::lsei(A=cbind(X0,u),B=Y0,E=e,F=f,G=g,H=h,type=2)$X
      w.hat <- c(temp[1:J])
      u.hat <- Y0-X0%*%w.hat
      sigma <- c(temp[J+1])
      return(list(u.hat=u.hat,w.hat=w.hat,sigma = sigma))
    }
  }
}



DRoSC <- function(Y0,Y1,X0,X1,lambda = NULL, intercept = F, nu = 0.01,eta = 0.01,
                 c = 1/2, 
                 Inference = F,alpha = 0.05, c.sample = 0.01,M = 500,
                 dependent = F,
                 true_mu = NULL,true_beta = NULL,
                 gamma = 1,alpha0 = 0.01){
  # method <- match.arg(method)
  T0 <- length(Y0)
  T1 <- length(Y1)
  N <- ncol(X0)
  Y0 <- as.vector(Y0)
  Y1 <- as.vector(Y1)
  
  if (length(lambda)==0) {
    lambda <- 0
  }
  
  # Normalize columns of X
  SC.fit <- sc(Y0,X0,intercept = intercept)
  w.hat <- SC.fit$w.hat
  sig <- as.numeric(sqrt(var(SC.fit$u.hat)))
  
  
  
  ## Pre-treatment period ##
  ## Uncertainty class ##
  X0.cov <- t(X0)%*%X0/T0 # Sigmahat
  # X0.norm <- max(min(eigen(X0.cov)$values),0.01)
  Y0X0.vec <- (Y0*X0)
  Y0X0.mean <- apply(Y0X0.vec,2,mean) # gammahat
  
  
  
  ## Post-treatment period
  ## Objective function ###
  EY1 <- mean(Y1) # muhat_Y
  EX1 <- apply(X1,2,mean) # muhat
  
  VarY1 <- ifelse(dependent,
                  as.numeric(HAC_cov_matrix(Y1,lag = floor(4 * (T1 / 100)^(2 / 9)))),
                  var(Y1))
  if (dependent) {
    CovX1 <-HAC_cov_matrix(X1, lag = floor(4 * (T1 / 100)^(2 / 9)))
  } else{
    CovX1 <-cov(X1)
  }
  
  # Tuning parameter #
  tauHat <- NA
  maxnorm <- max(apply(X0,2,function(x){sqrt(mean(x^2))}))
  # D <- diag(1/maxnorm)
  ### Sum-to-one ###
  e <- cbind(matrix(1,1,N))
  f <- 1
  check <- 0
  if (nu>0) {
    
    while(is.na(tauHat)){
      c0 <- log(max(T0,N))^c/sqrt(T0)*(nu*sig*maxnorm+eta*lambda)
      thres <- lambda+ c0
      # lambda.alpha <- lambda+c0
      g <- rbind(diag(x=1,N,N),X0.cov,-X0.cov)
      h <- rbind(matrix(0,N,1),cbind(c(Y0X0.mean-thres,-(Y0X0.mean+thres))))
      tryCatch(betaHat <-limSolve::lsei(A=EX1,B=EY1,E=e,F=f,G=g,H=h,type=2)$X,
               error = function(e) { betaHat <<- NA})
      tryCatch(tauHat<- EY1-EX1%*%betaHat,
               error = function(e) { tauHat <<-NA})
      nu <-  1.25*nu
      check <- check+1
    }
  } else{
    thres <- lambda+log(max(T0,N))^c/sqrt(T0)*(nu*sig*maxnorm+eta*lambda)
    # lambda.alpha <- lambda+c0
    g <- rbind(diag(x=1,N,N),X0.cov,-X0.cov)
    h <- rbind(matrix(0,N,1),cbind(c(Y0X0.mean-thres,-(Y0X0.mean+thres))))
    tryCatch(betaHat <-limSolve::lsei(A=EX1,B=EY1,E=e,F=f,G=g,H=h,type=2)$X,
             error = function(e) { betaHat <<- NA})
    tryCatch(tauHat<- EY1-EX1%*%betaHat,
             error = function(e) { tauHat <<-NA})
  }
  if (!is.na(tauHat)) {
    tauHat <- as.numeric(tauHat)
  }
  
  
  
  
  
  
  
  ###################
  #### Inference ####
  ###################
  if (Inference) {
    X0.cov.vec <- vecl(X0.cov) # vecl(\Sigma)
    vecl.X0t <- matrix(0,nrow = T0,ncol = length(X0.cov.vec))
    for (t in 1:T0) {
      vecl.X0t[t,] <- vecl(t(X0[t,,drop = F])%*%X0[t,,drop = F])
    }
    if (dependent) {
      cov.vecl.X0 <- HAC_cov_matrix(vecl.X0t, lag = floor(4 * (T0 / 100)^(2 / 9)))
    } else{
      cov.vecl.X0 <- cov(vecl.X0t)
    }
    
    
    dW <- max(gamma*max(abs(cov.vecl.X0)),1) 
    if (T0>N*(N+1)/2) {
      dW <- 0
    }
    Sigma.cov <- cov.vecl.X0/T0+dW/T0*diag(ncol(cov.vecl.X0))
    if (dependent) {
      Y0X0.cov <- HAC_cov_matrix(Y0X0.vec, lag = floor(4 * (T0 / 100)^(2 / 9)))
    } else{
      Y0X0.cov <- cov(Y0X0.vec)
    }
    
    
    dZ <- max(gamma*max(abs(Y0X0.cov)),1)
    if (T0>N) {
      dZ <- 0 
    }
    
    
    dX <- max(gamma*max(abs(CovX1)),1)
    if (T1>N) {
      dX <- 0
    }
    prop <- 0
    
    while(prop< 0.1){
      beta.mat <- matrix(NA,nrow = M,ncol = N)
      mu.mat <- matrix(NA,nrow = M,ncol = N)
      sample.XX <- list()
      sample.XY <- list()
      sample.X <- sample.Y <- list()
      for (m in 1:M) {
        skip_to_next <<- FALSE
        if (!cond) { # unconditional
          XX.repro <- MASS::mvrnorm(mu = X0.cov.vec,Sigma = Sigma.cov)
          max(abs(XX.repro-X0.cov.vec)/sqrt(diag(Sigma.cov)))
          # drop outlier
          if (max(abs(XX.repro-X0.cov.vec)/sqrt(diag(Sigma.cov)))>1.05*qnorm(1-alpha0/(2*length(XX.repro)))) {
            next
          }
          tmp <- matrix(0,nrow = N,ncol = N)
          tmp[lower.tri(tmp,diag = T)] <- XX.repro
          # Symmetric
          for(l in 2:N) {
            for(k in 1:(l-1)) {
              tmp[k,l] = tmp[l,k]
            }
          }
          Diag.XX<-diag(eigen(tmp)$values)
          for(ind in 1:N){
            Diag.XX[ind,ind]<-max(Diag.XX[ind,ind],1e-3)
          }
          ev <- eigen(tmp)
          XX.positive<-ev$vectors%*%Diag.XX%*%t(ev$vectors)
        } else{
          XX.positive <- X0.cov
        }
        sample.XX[[m]] <- XX.positive
        
        XY.repro <- MASS::mvrnorm(mu = Y0X0.mean,
                                  Sigma = Y0X0.cov/T0+dZ/T0*diag(ncol(Y0X0.cov)))
        if (max(abs(XY.repro-Y0X0.mean)/sqrt(diag(Y0X0.cov/T0+dZ/T0*diag(ncol(Y0X0.cov)))))>1.05*qnorm(1-alpha0/(2*length(XY.repro)))) {
          next
        }
        sample.XY[[m]] <- XY.repro
        EY1.m <- rnorm(1,mean = EY1, sd = sqrt((VarY1)/T1))
        if (max(abs(EY1.m-EY1)/sqrt((VarY1)/T1))>1.05*qnorm(1-alpha0/(2*length(EY1.m)))) {
          next
        }
        if (cond) {
          EX1.m <- EX1
        } else{
          EX1.m <- MASS::mvrnorm(n=1,mu = EX1,Sigma = CovX1/T1+dX/T1*diag(N))
          if (max(abs(EX1.m-EX1)/sqrt(diag(CovX1/T1+dX/T1*diag(N))))>1.05*qnorm(1-alpha0/(2*length(EX1.m)))) {
            next
          }
        }
        sample.X[[m]] <- EX1.m
        sample.Y[[m]] <- EY1.m
        A <- (EX1.m)%*%t(EX1.m)
        beta.m <- NA
        lambda.hat <- lambda+c.sample/sqrt(T0)*(log(min(T0,T1))/M)^{1/(1+(N*(N+5))/2)}*max(sqrt(diag(Y0X0.cov/T0+dZ/T0*diag(ncol(Y0X0.cov)))),
                                                                                  sqrt(diag(Sigma.cov)))

        g <- rbind(diag(x=1,N,N),XX.positive,-XX.positive)
        h <- rbind(matrix(0,N,1),cbind(c(XY.repro-lambda.hat,-(XY.repro+lambda.hat))))
        tryCatch(beta.m <-limSolve::lsei(A=EX1.m,B=EY1.m,E=e,F=f,G=g,H=h,type=2)$X,
                 error = function(e) { skip_to_next <<- TRUE})
        beta.mat[m,] <- beta.m
        mu.mat[m,] <- EX1.m
      }
      valid_rows <- rowSums(!is.na(beta.mat)) > 0
      prop <- mean(valid_rows)
      if (prop<0.1) {
        c.sample <- 1.25*c.sample
      }
    }
    mu.valid <- mu.mat[valid_rows,]
    beta.valid <- beta.mat[valid_rows,]
    sample.XX <- sample.XX[valid_rows]
    sample.XY <-sample.XY[valid_rows]
    sample.X <- sample.X[valid_rows]
    sample.Y <- sample.Y[valid_rows]
    tau.vec <- vector()
    Int.mat <- matrix(NA,nrow = nrow(beta.valid),ncol = 2)
    SE <- sqrt(VarY1)
    for (i in 1:nrow(beta.valid)) {
      tau.m <- EY1 - EX1.m%*%beta.valid[i,]
      tau.vec[i] <- tau.m
      # SE <- sd(Y1-EY1)
      Int.mat[i,] <- c(tau.m-qnorm(1-(alpha-alpha0)/2)*SE/sqrt(T1),tau.m+qnorm(1-(alpha-alpha0)/2)*SE/sqrt(T1))
    }
    if (!is.null(true_mu) & !is.null(true_beta)) {
      mu.diff = sweep(mu.valid, 2, true_mu, FUN = "-")
      beta.diff = sweep(beta.valid, 2, true_beta, FUN = "-")
      negli <- abs(-mu.diff%*%true_beta) + abs(rowSums(-mu.valid*beta.diff))
      
      # Find the row index with the minimum mu1\tr(beta.m-beta.star)
      min_index <- which.min(negli)
      
      beta.min <- beta.valid[min_index,]
      mu.min <- mu.valid[min_index,]
      negli.min <- min(negli)
    }
    if (all(is.na(Int.mat))) {
      CI.tau = NA
    } else{
      suppressWarnings({
        uni = intervals::Intervals(Int.mat)
        CI.tau = as.matrix(intervals::interval_union(uni))
      })
    }
  }
  
  if (Inference) {
    if (!is.null(true_beta) & !is.null(true_mu)) {
      return(list(betaHat = betaHat,tauHat = tauHat,check = check,
                  CI.tau=CI.tau,Int.mat = Int.mat,
                  mu.min = mu.min,beta.min = beta.min, negli.min =  negli.min,
                  prop = prop,c.sample = c.sample ))
    } else{
      return(list(betaHat = betaHat,tauHat = tauHat,check = check,
                  CI.tau=CI.tau,prop = prop,c.sample = c.sample))
    }
    
  } else{
    return(list(betaHat = betaHat,tauHat = tauHat,check = check))
  }
  
  
}