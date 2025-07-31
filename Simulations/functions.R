


### anciliary functions
fdp_power <- function(selected_index, signal_index){
  num_selected <- length(selected_index)
  tp <- length(intersect(selected_index, signal_index))
  fp <- num_selected - tp
  fdp <- fp / max(num_selected, 1)
  power <- tp / length(signal_index)
  return(list(fdp = fdp, power = power))
}


Sigmoid = function(x){
  1/(1+exp(-x))
}



analysis <- function(mm, ww, q){
    cutoff_set = max(ww)
    for(t in ww){
        ps = length(mm[mm>=t])
        ng = length(na.omit(mm[mm<=-t]))
        rto = (ng)/max(ps, 1)
        if(rto<=q){
            cutoff_set = c(cutoff_set, t)
        }
    }
    cutoff = min(cutoff_set)
    selected_index = which(mm >= cutoff)

    return(selected_index)
}


AdaFDR = function(x, y, family = "gaussian", off = 1, reduce = 1, lambda_min = 0,
                    kappa = seq(1,10,1), sigma = 1, fdr = 0.1){
  ## Input parameters
  # x: the n by p design matrix
  # y: the n by 1 response vector
  # family: 'gaussian' or 'binomial'
  # off: control the bias when approximing FDP, if off = 0, the bias = 0.
  # reduce: rheduce = 1 denotes get the feature importance scores 
  #   on the retained set of initial LASSO estimate
  # kappa: the parameter in weighting vector, kappa>=1
  # sigma: the noise level by Gaussian perturbation
  # fdr: the pre-designed FDR level

  ## Out put
  # selected: the selected variable set
  # Statistics: the statistics associated with selected variable set
  # kappa_max: the optimal parameter in weighting vector

  n = dim(x)[1] # sample size
  p = dim(x)[2] # dimension

  ## fit the initial lasso regression on (x, y)
  fit0.lasso <- cv.glmnet(x = x, y = y, family = family, alpha = 1, standardize=F)
  if(lambda_min == 1){
    # estimated tuning parameter
    lambda_est <- fit0.lasso$lambda.min
    # estimated regression coefficient
    initial_beta <- as.vector(coef(fit0.lasso, s = 'lambda.min'))[-1]
  }else{
    # estimated tuning parameter
    lambda_est <- fit0.lasso$lambda.1se
    # estimated regression coefficient
    initial_beta <- as.vector(coef(fit0.lasso, s = 'lambda.1se'))[-1]
  }
  if(reduce == 1){
      # selected variable set
      S1 <- which(initial_beta != 0)
      p_hat = length(S1)
  }else{
      # selected variable set
      S1 <- c(1:p)
      p_hat = length(S1)
  }
  # the feature importance score vector
  M_hat <- M_tilde <- numeric(p_hat)
  X_hat = X[ ,S1]
  for(j_var in 1:p_hat){
    # generate gaussian noise
    noise1 <- rnorm(n, 0, sigma)
    noise2 <- rnorm(n, 0, sigma)
    # combine new data
    Xnew1 <- cbind(X_hat[ ,j_var] + noise1, 
                    X_hat[ ,j_var] - noise1, X_hat[ ,-j_var])
    Xnew2 <- cbind(X_hat[ ,j_var] + noise2, 
                    X_hat[ ,j_var] - noise2, X_hat[ ,-j_var])
    # fit lasso regression
    fit1.lasso <- glmnet(x = Xnew1, y = y, family = family, 
                          lambda = lambda_est, alpha = 1)    
    fit2.lasso <- glmnet(x = Xnew2, y = y, family = family, 
                          lambda = lambda_est, alpha = 1)
    # get estimated regression coefficient
    beta_hat1 <- as.vector(coef(fit1.lasso))[-1]
    beta_hat2 <- as.vector(coef(fit2.lasso))[-1]
    
    b_hat <- abs(2*beta_hat1[1:2])
    b_tilde <- abs(2*beta_hat2[1:2])
    # get variable importance score
    M_hat[j_var] <- (b_hat[1]*b_hat[2])
    M_tilde[j_var] <- (b_tilde[1]*b_tilde[2])
    
  }
  
  # the length of kappa
  len_kappa = length(kappa)
  # the weight vector for different kappa
  weight = matrix(1, p, len_kappa)
  for(k in 1:len_kappa){
    weight[S1, k] = 1/(abs(initial_beta[S1]) + 1 )^(kappa[k])
  }
  # get test statistics
  M_pgm <- matrix(0, p, len_kappa)
  for(k in 1:len_kappa){
    M_pgm[S1, k] = M_hat - weight[S1, k]*M_tilde
  }
  
  # the index of selected variables
  selected_index = list()
  # the number of selected variables
  len_selected = numeric(len_kappa)

  ## find the threshold
  for(k in 1:len_kappa){
    M_pgm_k = M_pgm[,k]
    ww = abs(M_pgm_k)
    cutoff_set <- max(ww)
    bias = numeric(p)
    for(j in 1:p){
      t <- ww[j]
      # the selected variables given t
      ps <- length(M_pgm_k[ M_pgm_k > t ])
      # the false selected variables given t
      ng <- length(na.omit(M_pgm_k[M_pgm_k < -t]))
      ind1 = which(M_pgm_k < t/weight[, k])
      ind2 = which(M_pgm_k < t)
      # the bias term in the estimated FDP
      bias[j] = length(ind1) - length(ind2)
      if(off == 1){
        rto <- (ng + bias[j])/ max(ps, 1)
      }else{
       rto <- (ng)/ max(ps, 1)
      }
      if(rto <= FDR){
        cutoff_set <- c(cutoff_set, t)
      }
    }
    # threshold value
    cutoff <- min(cutoff_set) 
    # the selected variable set given kappa
    selected_index[[k]] <- which( M_pgm_k > cutoff)
    # the number of selected variables given kappa
    len_selected[k] = length(selected_index[[k]])
 }
 # return the index of optimal kappa 
 k_max = which.max(len_selected)
 # the selected variables associated with optimal kappa
 selected_index = selected_index[[k_max]]
 M_pgm_max = M_pgm[, k_max]
 kappa_max = kappa[k_max]
 return(list(selected = selected_index, Statistics = M_pgm_max, 
        kappa_max = kappa_max))
}


####### LM data-splitting methods (DS and MDS)

f = function(x){
  exp(x)/(1+exp(x))
}

#Calculate the precision matrix
get_M = function(data){
  ## data: design matrix
  n = dim(data)[1]
  p = dim(data)[2]
  C = matrix(0, nrow = p, ncol = p)
  tau_square = numeric(p)
  ## Nodewise Lasso regression
  for(j in 1:p){
    y = data[,j]
    X <- data[, -j]
    ## To save computation time, we only do the cv once
    if(j == 1){
      cvfit <- cv.glmnet(X, y, nfolds = 10, nlambda = 200, intercept = F, standardize = F)
      ## We will use this same lambda for the following Lasso regression
      lambda <- cvfit$lambda.min
    }
    beta1 <- as.vector(glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = F, standardize = F)$beta)
    C[j, -j] = -beta1
    tau_square[j] = mean((y-X%*%beta1)*y)
  }
  diag(C) = 1
  T_mat = diag(tau_square)
  M = solve(T_mat)%*%C
  M
}


DS_lm = function(x, y,q){
  ## x: design matrix
  ## y: response variable
  n = dim(x)[1]; p = dim(x)[2]
  ## split the data into two halves and run Lasso
  sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
  sample_index2 <- setdiff(c(1:n), sample_index1)
  X1 = x[sample_index1, ]
  y1 = y[sample_index1]
  X2 = x[sample_index2, ]
  y2 = y[sample_index2]
  
  cvfit <- cv.glmnet(x, y, family = "gaussian", type.measure = "mse", nfolds = 10, nlambda = 200,  intercept = F, standardize = F)
  lambda <- cvfit$lambda.min/sqrt(2)
  beta1 <- as.vector(glmnet(X1, y1, family = "gaussian", alpha = 1, lambda = lambda, intercept = F, standardize = F)$beta)
  beta2 <- as.vector(glmnet(X2, y2, family = "gaussian", alpha = 1, lambda = lambda, intercept = F, standardize = F)$beta)
  ## Calculate the weight for each observation
  w1 = as.vector(sqrt(f(X1%*%beta1)*(1-f(X1%*%beta1))))
  w2 = as.vector(sqrt(f(X2%*%beta2)*(1-f(X2%*%beta2))))
  ## Get the precision matrix for the reweighted design matrix
  X1_beta = diag(w1)%*%X1
  X2_beta = diag(w2)%*%X2
  M1 = get_M(X1_beta)
  M2 = get_M(X2_beta)
  ## Get the debiased Lasso estimator
  beta1_d = as.numeric(beta1 + 2/n*M1%*%t(x[sample_index1, ])%*%(y[sample_index1] - f(x[sample_index1,] %*% beta1)))
  beta2_d = as.numeric(beta2 + 2/n*M2%*%t(x[sample_index2, ])%*%(y[sample_index2] - f(x[sample_index2,] %*% beta2)))
  ## Get the variance for the debiased Lasso estimator (up to a scaling factor)
  sigma1 = sqrt(diag(M1%*%(t(X1)%*%diag(diag((y1-f(X1%*%beta1))%*%t(y1-f(X1%*%beta1))))%*%X1/(n/2))%*%t(M1)))
  sigma2 = sqrt(diag(M2%*%(t(X2)%*%diag(diag((y2-f(X2%*%beta2))%*%t(y2-f(X2%*%beta2))))%*%X2/(n/2))%*%t(M2)))
  M = beta1_d*beta2_d/(sigma1*sigma2)
  ## Get the selected_index
  select_index = analysis(M, abs(M), q)
  num_select <- length(select_index)
  return(list(select_index = select_index, num_select = num_select))
}


MDS_lm = function(x, y, num_split, q){
  ## x: design matrix
  ## y: response variable
  ## num_split: The number of split
  n = dim(x)[1]; p = dim(x)[2]
  inclusion_rate_multiple <- matrix(0, nrow = num_split, ncol = p)
  fdr_multiple <- rep(0, num_split)
  power_multiple <- rep(0, num_split)
  num_select <- rep(0, num_split)
  ## Run DS num_split times.
  for(iter in 1:num_split){
    ## The same code as DS
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)
    X1 = x[sample_index1, ]
    y1 = y[sample_index1]
    X2 = x[sample_index2, ]
    y2 = y[sample_index2]
    
    cvfit <- cv.glmnet(x, y, family = "gaussian", type.measure = "mse", nfolds = 10, nlambda = 200, intercept = F, standardize = F)
    lambda <- cvfit$lambda.min/sqrt(2)
    beta1 <- as.vector(glmnet(X1, y1, family = "gaussian", alpha = 1, lambda = lambda, intercept = F, standardize = F)$beta)
    beta2 <- as.vector(glmnet(X2, y2, family = "gaussian", alpha = 1, lambda = lambda, intercept = F, standardize = F)$beta)
    
    w1 = as.vector(sqrt(f(X1%*%beta1)*(1-f(X1%*%beta1))))
    w2 = as.vector(sqrt(f(X2%*%beta2)*(1-f(X2%*%beta2))))
    X1_beta = diag(w1)%*%X1
    X2_beta = diag(w2)%*%X2
    M1 = get_M(X1_beta)
    M2 = get_M(X2_beta)
    beta1_d = as.numeric(beta1 + 2/n*M1%*%t(x[sample_index1, ])%*%(y[sample_index1] - f(x[sample_index1,] %*% beta1)))
    beta2_d = as.numeric(beta2 + 2/n*M2%*%t(x[sample_index2, ])%*%(y[sample_index2] - f(x[sample_index2,] %*% beta2)))
    sigma1 = sqrt(diag(M1%*%(t(X1)%*%diag(diag((y1-f(X1%*%beta1))%*%t(y1-f(X1%*%beta1))))%*%X1/(n/2))%*%t(M1)))
    sigma2 = sqrt(diag(M2%*%(t(X2)%*%diag(diag((y2-f(X2%*%beta2))%*%t(y2-f(X2%*%beta2))))%*%X2/(n/2))%*%t(M2)))
    M = beta1_d*beta2_d/(sigma1*sigma2)
    #Get the selected index for a single split
    current_selected_index = analysis(M, abs(M), q)
    num_select[iter] <- length(current_selected_index)
    inclusion_rate_multiple[iter, current_selected_index] <- 1/num_select[iter]
    result = fdp_power(current_selected_index, signal_index)
    fdr_multiple[iter] = result$fdp
    power_multiple[iter] = result$power
  }
  ## single splitting result
  single_split_fdr <- fdr_multiple[1]
  single_split_power <- power_multiple[1]
  ## multiple splitting result
  inclusion_rate <- apply(inclusion_rate_multiple, 2, mean)
  ## sort the features according to the inclusion rate
  feature_rank <- order(inclusion_rate)
  ## discard the features with zero inclusion rate
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  null_variable <- numeric()
  ## Choose cutoff
  for(feature_index in 1:length(feature_rank)){
    if(sum(inclusion_rate[feature_rank[1:feature_index]]) > q){
      break
    }else{
      null_variable <- c(null_variable, feature_rank[feature_index])
    }
  }
  select_index <- setdiff(feature_rank, null_variable)
  num_select <- length(select_index)
  return(list(select_index = select_index, num_select = num_select))
}



### GLM data-splitting methods (DS and MDS)
DS_MDS_glm <- function(X, y, signal_index, num_split = 50,  q = 0.1){
  n <- dim(X)[1]; p <- dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  fdp <- rep(0, num_split)
  power <- rep(0, num_split)
  num_select <- rep(0, num_split)
  
  for(iter in 1:num_split){
    ### randomly split the data
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)
    ### get the penalty lambda for Lasso
    fit1 <- glm(y[sample_index1] ~ X[sample_index1,] - 1, family = 'binomial')
    beta1 <- fit1$coefficients
    nonzero_index <- which(beta1 != 0)
    if(length(nonzero_index)!=0){
      ### run OLS on the second half of the data, restricted on the selected features
      beta2 <- rep(0, p)
      fit2 <- glm(y[sample_index2] ~ X[sample_index2,] - 1, family = 'binomial')
      beta2 <- fit2$coefficients
      ## Calculate the variance of each estimator(up to a scaling factor)
      tau1 <- 1/sqrt(diag(ginv(t(X[sample_index1, ]) %*% X[sample_index1, ]))) #solve
      tau2 <- 1/sqrt(diag(ginv(t(X[sample_index2, ]) %*% X[sample_index2, ])))
      M <- beta1*beta2*tau1*tau2
      selected_index <- analysis(M, abs(M), q)
      
      ### number of selected variables
      if(length(selected_index)!=0){
        num_select[iter] <- length(selected_index)
        inclusion_rate[iter, selected_index] <- 1/num_select[iter]
        
        ### calculate fdp and power
        result <- fdp_power(selected_index, signal_index)
        fdp[iter] <- result$fdp
        power[iter] <- result$power
      }
    }
  }
  ### single data-splitting (DS) result
  DS_fdp <- fdp[1]
  DS_power <- power[1]
  
  ### multiple data-splitting (MDS) result
  inclusion_rate <- apply(inclusion_rate, 2, mean)
  ### rank the features by the empirical inclusion rate
  feature_rank <- order(inclusion_rate)
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  if(length(feature_rank)!=0){
    null_feature <- numeric()
    ### backtracking
    for(feature_index in 1:length(feature_rank)){
      if(sum(inclusion_rate[feature_rank[1:feature_index]]) > q){
        break
      }else{
        null_feature <- c(null_feature, feature_rank[feature_index])
      }
    }
    selected_index <- setdiff(feature_rank, null_feature)
    ### calculate fdp and power
    result <- fdp_power(selected_index, signal_index)
    MDS_fdp <- result$fdp
    MDS_power <- result$power
  }
  else{
    MDS_fdp <- 0
    MDS_power <- 0
  }
  return(list(DS_fdp = DS_fdp, DS_power = DS_power, MDS_fdp = MDS_fdp, MDS_power = MDS_power))
}



DS_glm = function(x, y, signal_index, q){
  ## x: design matrix
  ## y: response variable
  n = dim(x)[1]; p = dim(x)[2]
  ## split the data into two halves and run Lasso
  sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
  sample_index2 <- setdiff(c(1:n), sample_index1)
  X1 = x[sample_index1, ]
  y1 = y[sample_index1]
  X2 = x[sample_index2, ]
  y2 = y[sample_index2]
  
  cvfit <- cv.glmnet(X1, y1, family = "poisson", type.measure = "mse", nfolds = 10, intercept = F, standardize = F)
  lambda <- cvfit$lambda.min
  beta1 = as.vector(coef(cvfit, lambda))[-1]
  cvfit <- cv.glmnet(X2, y2, family = "poisson", type.measure = "mse", nfolds = 10, intercept = F, standardize = F)
  lambda <- cvfit$lambda.min
  beta2 = as.vector(coef(cvfit, lambda))[-1]
  
  ## Calculate the weight for each observation
  w1 = as.vector(sqrt(f(X1%*%beta1)))
  w2 = as.vector(sqrt(f(X2%*%beta2)))
  
  ## Get the precision matrix for the reweighted design matrix
  X1_beta = diag(w1)%*%X1
  X2_beta = diag(w2)%*%X2
  M1 = get_M(X1_beta)
  M2 = get_M(X2_beta)
  
  ## Get the debiased Lasso estimator
  beta1_d = as.numeric(beta1 + 2/n*M1%*%t(x[sample_index1, ])%*%(y[sample_index1] - f(x[sample_index1,] %*% beta1)))
  beta2_d = as.numeric(beta2 + 2/n*M2%*%t(x[sample_index2, ])%*%(y[sample_index2] - f(x[sample_index2,] %*% beta2)))
  
  ## Get the variance for the debiased Lasso estimator (up to a scaling factor)
  sigma1 = sqrt(diag(M1%*%(t(X1)%*%diag(diag((y1-f(X1%*%beta1))%*%t(y1-f(X1%*%beta1))))%*%X1/(n/2))%*%t(M1)))
  sigma2 = sqrt(diag(M2%*%(t(X2)%*%diag(diag((y2-f(X2%*%beta2))%*%t(y2-f(X2%*%beta2))))%*%X2/(n/2))%*%t(M2)))
  M = beta1_d*beta2_d/(sigma1*sigma2)
  
  ## Get the selected_index
  selected_index = analysis(M, abs(M), q)
  
  ## Calculate fdp and power
  result = fdp_power(selected_index, signal_index)
  DS_fdp <- result$fdp
  DS_power <- result$power
  return(list(DS_fdp = DS_fdp, DS_power = DS_power))
}



MDS_glm = function(x, y, signal_index, num_split, q){
  ## x: design matrix
  ## y: response variable
  ## num_split: The number of split
  n = dim(x)[1]; p = dim(x)[2]
  inclusion_rate_multiple <- matrix(0, nrow = num_split, ncol = p)
  fdr_multiple <- rep(0, num_split)
  power_multiple <- rep(0, num_split)
  num_select <- rep(0, num_split)
  ## Run DS num_split times.
  for(iter in 1:num_split){
    ## The same code as DS
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)
    X1 = x[sample_index1, ]
    y1 = y[sample_index1]
    X2 = x[sample_index2, ]
    y2 = y[sample_index2]
    
    cvfit <- cv.glmnet(X1, y1, family = "poisson", type.measure = "mse", nfolds = 10, intercept = F, standardize = F)
    lambda <- cvfit$lambda.min
    beta1 = as.vector(coef(cvfit, lambda))[-1]
    cvfit <- cv.glmnet(X2, y2, family = "poisson", type.measure = "mse", nfolds = 10, intercept = F, standardize = F)
    lambda <- cvfit$lambda.min
    beta2 = as.vector(coef(cvfit, lambda))[-1]
    
    w1 = as.vector(sqrt(f(X1%*%beta1)))
    w2 = as.vector(sqrt(f(X2%*%beta2)))
    
    X1_beta = diag(w1)%*%X1
    X2_beta = diag(w2)%*%X2
    M1 = get_M(X1_beta)
    M2 = get_M(X2_beta)
    beta1_d = as.numeric(beta1 + 2/n*M1%*%t(x[sample_index1, ])%*%(y[sample_index1] - f(x[sample_index1,] %*% beta1)))
    beta2_d = as.numeric(beta2 + 2/n*M2%*%t(x[sample_index2, ])%*%(y[sample_index2] - f(x[sample_index2,] %*% beta2)))
    sigma1 = sqrt(diag(M1%*%(t(X1)%*%diag(diag((y1-f(X1%*%beta1))%*%t(y1-f(X1%*%beta1))))%*%X1/(n/2))%*%t(M1)))
    sigma2 = sqrt(diag(M2%*%(t(X2)%*%diag(diag((y2-f(X2%*%beta2))%*%t(y2-f(X2%*%beta2))))%*%X2/(n/2))%*%t(M2)))
    M = beta1_d*beta2_d/(sigma1*sigma2)
    #Get the selected index for a single split
    current_selected_index = analysis(M, abs(M), q)
    num_select[iter] <- length(current_selected_index)
    inclusion_rate_multiple[iter, current_selected_index] <- 1/num_select[iter]
    result = fdp_power(current_selected_index, signal_index)
    fdr_multiple[iter] = result$fdp
    power_multiple[iter] = result$power
  }
  ## single splitting result
  single_split_fdr <- fdr_multiple[1]
  single_split_power <- power_multiple[1]
  ## multiple splitting result
  inclusion_rate <- apply(inclusion_rate_multiple, 2, mean)
  ## sort the features according to the inclusion rate
  feature_rank <- order(inclusion_rate)
  ## discard the features with zero inclusion rate
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  null_variable <- numeric()
  ## Choose cutoff
  for(feature_index in 1:length(feature_rank)){
    if(sum(inclusion_rate[feature_rank[1:feature_index]]) > q){
      break
    }else{
      null_variable <- c(null_variable, feature_rank[feature_index])
    }
  }
  selected_index <- setdiff(feature_rank, null_variable)
  ## Calculate fdp and power
  result = fdp_power(selected_index, signal_index)
  MDS_fdp <- result$fdp
  MDS_power <- result$power
  return(list(MDS_fdp = MDS_fdp, MDS_power = MDS_power))
  
}





