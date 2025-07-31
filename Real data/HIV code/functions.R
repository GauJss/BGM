###### design matrix function
design_matrix <- function(data, col, row){
  dataX=data[,col:ncol(data)]
  dataX=dataX[row:nrow(dataX),]
  for (i in 1:ncol(dataX)) {
    dataX[which(dataX[,i]=="-"),i]=0
    dataX[which(dataX[,i]=="."),i]=0
    dataX[which(dataX[,i]==""),i]=0
  }
  design_X=matrix(0,nrow = nrow(dataX),ncol = 1)
  for (i in 1:ncol(dataX)) {
    xxdata=dataX[,i]
    xrow=rownames(table(xxdata))
    xrow=xrow[which(xrow!="0")] 
    ll=length(xrow)
    if (ll>0){
      pix=matrix(rep(xxdata,length(xrow)),ncol=length(xrow))
      for (j in 1:ncol(pix)) {
        tt=which(pix[,j]==xrow[j])
        pix[tt,j]=1
        pix[tt,-j]=0
      }
      xrowname=paste(as.character(i),xrow,sep = "")
      colnames(pix)=xrowname
      design_X=cbind(design_X,pix)
    }}
  design_X=design_X[,-1]
  return(design_X)
}


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
        ps = length(mm[mm>=1])
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
  X_hat = x[ ,S1]
  X_hat=as.data.frame(X_hat)
  X_hat=sapply(X_hat,as.numeric)
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
      if(rto <= fdr){
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







### data-splitting methods (DS and MDS)
DS_glm <- function(X, y, signal_index, num_split = 50,  q = 0.1){
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

    cvfit <- cv.glmnet(X[sample_index1, ], y[sample_index1], family = "binomial")
    lambda <- cvfit$lambda.min
    
    ### run Lasso on the first half of the data
    beta1 <- as.vector(glmnet(X[sample_index1, ], y[sample_index1], family = "binomial", alpha = 1, lambda = lambda)$beta)
    nonzero_index <- which(beta1 != 0)
    if(length(nonzero_index)!=0){
      ### run OLS on the second half of the data, restricted on the selected features
      beta2 <- rep(0, p)
      beta2[nonzero_index] <- as.vector(glm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1, family = "binomial" )$coeff)
      
      ### calculate the mirror statistics
      #M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
      M <- abs(beta1 + beta2) - abs(beta1 - beta2)
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

