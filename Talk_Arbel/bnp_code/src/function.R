####################################################################################
# code by Riccardo Corradin
# Based on ongoing work:
# Arbel, J., Corradin, R., Nipoti, B. (2017). 
# Robust Bayesian nonparametric methods for density estimation 
# and clustering on the phase-space. In preparation.
####################################################################################

#####################
# Required packages #
#####################

library(mvtnorm)
library(abind)

#####################
# Multivariate code #
#####################

# data = multivariate dataset
# grid = grid of points where evaluate the density
# nsim = number of total iteration
# nburn = number of burn-in iteration

# theta = initial value for the mass parameter of DP
# m0 = initial value for the mean of location component of base measure
# k0 = initial value for correction parameter for the dispersion of location component
# nu0 = gdl of Inverse Wishart distribution of scale component of base measure
# Lambda0 = initial value for the matrix which characterize the Inverse Wishart distribution

# t1, t2 = hyperparameters for Gamma hyperdistribution on theta
# m1, m2 = hyperparameters (location - scale) for Normal hyperdistribution on m0
# b1, b2 = hyperparameters for Wishart hyperdistribution on Lambda0
# k1, k2 = hyperparameters for Gamma hyperdistribution on k0


MCMCcoolMulti <- function(data, grid, nsim, nburn, theta, m0, 
                          k0, nu0, Lambda0, t1, t2, m1, m2, b1, b2, k1, k2){
  
  # define n as number of observation
  # define d as number of variables
  as.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  ng <- nrow(grid)
  
  # initialize conf for save each configuration for single iteration
  # initialize theta for trace the value of the DP parameter
  conf <- matrix(0, nrow = nsim, ncol = n)
  theta_out <- c()
  
  # split randomly the observations into 2 groups 
  # initialize the param object with the parameters for each group
  #    NB: 1 extra dimension for mantain the object as tensor in the code
  # initialize the matrix to safe all iterations computed density
  #    on the points given by grid object
  S <- sample(c(1,2), size = n, replace = T)
  param <- abind(matrix(0, ncol = d + 1, nrow = d), matrix(cbind(rep(0, d), diag(rep(1,d))), ncol = d + 1), 
                 matrix(cbind(rep(0, d), diag(rep(1,d))), ncol = d + 1), along = 3)
  dest <- rep(0, nrow(grid))
  
  # iterate for the number of simulations
  pb <- txtProgressBar(0, nsim, style = 3)
  for(i in 1:nsim){
    
    # for each simulation
    # iterate over the observation
    for(j in 1:n){
      
      # select just the parameters needed 
      param <- param[,,c(1, sort(unique(S)) + 1)]
      S <- as.numeric(as.factor(S))
      
      # generate the probability vector 
      # loop from 2 over the param object
      # the index 1 is for the null value which preserve tenso structure
      # the probability for each block is proportional to the block size and the kernel
      # the probability for a new block is proportional to the prior predictive distribution
      # and  the mass parameter of the DP
      # if the block is new, add the parameters to for the new block in the object param
      p <- c()
      for(k in 2:dim(param)[3]){
        p[k - 1] <- sum(S[-j] == k - 1) * dmvnorm(x = data[j,], mean = param[,1,k], sigma = param[,-1,k], log = F)
        
      }
      p <- c(theta * dmvt(x = data[j,], df = nu0 - d + 1, delta = m0, sigma = Lambda0 * 
                            (k0 + 1) / (k0 *(nu0 -d + 1)), log = F), p)
      p <- p/sum(p)
      
      S[j] <- sample(x = 0:max(S), size = 1, prob = p)
      if(S[j] == 0){
        
        S[j] <- max(S) + 1
        param <- abind(param, cbind((k0 * m0 + data[j,]) / (k0 + 1), 
                                    (Lambda0 + (k0 / (k0 + 1) ) * (data[j,] - m0) %*% t(data[j,] - m0)) * 
                                      (k0 + 1) / (k0 *(nu0 -d + 1))) , along = 3)
        
      }
      
    }
    
    # discard useless parameters
    param <- param[,,c(1, sort(unique(S)) + 1)]
    S <- as.numeric(as.factor(S))
    conf[i,] <- S
    
    # for each object in param
    for(k in 2:dim(param)[3]){
      
      # update hyperparameters according to the proper formula
      # splitting the cases with 1 or more observation in the block
      if(sum(S == k - 1) > 1){
        
        kn <- k0 + sum(S == k - 1)
        nun <- nu0 + sum(S == k - 1)
        mn <- (k0 * m0 + colSums(data[S == k - 1,])) / kn
        Lambdan <- Lambda0 + (sum(S == k - 1) - 1) * var(data[S == k - 1,]) + 
          (k0 * sum(S == k - 1) / kn ) * (colMeans(data[S == k - 1,]) - m0) %*% t(colMeans(data[S == k - 1,]) - m0)
        
      }else{
        
        kn <- k0 + 1
        nun <- nu0 + 1
        mn <- (k0 * m0 + data[S == k - 1,]) / kn
        Lambdan <- Lambda0 + (k0 / kn ) * (data[S == k - 1,] - m0) %*% t(data[S == k - 1,] - m0)
        
      }
      
      # Update parameters from posterior distributions
      param[,-1,k] <- solve(rWishart(1, nun, solve(Lambdan))[,,1])
      param[,1,k] <- rmvnorm(1, mean = mn, sigma = param[,-1,k]/kn)
      
    }
    
    # update hyperparameters according to hyperprior distributions
    # split case for 1 or more blocks 
    # then compute the density for the current iteration and 
    # save it in the object dest
    if(dim(param)[3] > 2){
      
      k1n <- k1 + dim(param)[3] - 1
      k2n <- k2 + sum(apply(param[,,-1], 3, function(x) t(x[,1] - m0) %*% x[,-1] %*% (x[,1] - m0)))
      k0 <- rgamma(1, shape = k1n/2, rate = k2n/2)
      
      m2n <- solve(solve(m2) + k0 * matrix(rowSums(apply(param[,-1,-1],3,solve)), ncol = d, byrow = T))
      m1n <- m2n %*% (solve(m2) %*% m1 + k0 * rowSums(apply(param[,,-1],3, function(x) solve(x[,-1]) %*% x[,1])))
      m0 <- as.vector(rmvnorm(1, mean = m1n, sigma = m2n))
      
      b1n <- b1 + (dim(param)[3] - 1) * nu0
      b2n <- solve(solve(b2) + matrix(rowSums(apply(param[,-1,-1],3,solve)), ncol = d, byrow = T))
      Lambda0 <- rWishart(1, df = b1n, Sigma = b2n)[,,1]
      
      eta <- rbeta(1, theta + 1, n)
      pre <- (t1 + dim(param)[3] - 2) / ( t1 + dim(param)[3] - 2 + n * (t2 - log(eta)))
      theta <- rgamma(1, shape = t1 + dim(param)[3] - 1 - sample(0:1, 1, prob = c(1 - pre, pre)), rate = t2 - log(eta))
      
      if(i > nburn){
        dest <- dest + rowSums(apply(param[,,-1],3, function(x) dmvnorm(grid, mean = x[,1], sigma = x[,-1], log = F)) * 
                                 (matrix(rep(table(S)/n, ng), nrow = ng, byrow = T)))
      }
      
    }else{
      
      k1n <- k1 + 1
      k2n <- k2 + t(param[,1,-1] - m0) %*% param[,-1,-1] %*% (param[,1,-1] - m0)
      k0 <- rgamma(1, shape = k1n/2, rate = k2n/2)
      
      m2n <- solve(solve(m2) + k0 * solve(param[,-1,-1]))
      m1n <- m2n %*% (solve(m2) %*% m1 + k0 * solve(param[,-1,-1]) %*% param[,1,-1])
      m0 <- as.vector(rmvnorm(1, mean = m1n, sigma = m2n))
      
      b1n <- b1 + nu0
      b2n <- solve(solve(b2) + solve(param[,-1,-1]), ncol = d, byrow = T)
      Lambda0 <- rWishart(1, df = b1n, Sigma = b2n)[,,1]
      
      eta <- rbeta(1, theta + 1, n)
      pre <- t1 / ( t1 + n * (t2 - log(eta)))
      theta <- rgamma(1, shape = t1 + 1 - sample(0:1, 1, prob = c(1 - pre, pre)), rate = t2 - log(eta))
      
      if(i > nburn){
        dest <- dest + dmvnorm(grid, mean = param[,1,-1], sigma = param[,-1,-1], log = F) 
      }
      
    }
    # save the mass parameter of DP
    theta_out[i] <- theta         
    setTxtProgressBar(pb, i)
    
  }
  
  # return a list of 3 objects:
  # -the density object
  # -the configurations of the blocks for each iteration 
  # -the mass parameter object
  close(pb)
  return(list(dest/(nsim - nburn), conf[-c(1:nburn),], theta_out[-c(1:nburn)]))
}

##########################
# Bad 2d integral approx #
##########################

really_bad_int_approx <- function(grid, output){
  xval <- unique(grid[,1])
  yval <- unique(grid[,2])
  zval <- matrix(output[[1]], ncol = 40)
  
  for(j in 2:(length(xval)-1)){
    for(k in 2:length(yval)){
      area <- abs((xval[j] + xval[j-1])/2 - (xval[j+1] + xval[j])/2)  * abs((yval[j] + yval[j-1])/2 - (yval[j+1] + yval[j])/2) 
      inte <- inte + area * zval[j,k]
    }
  }
  return(inte)
}