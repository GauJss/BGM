
### simulation results with logistic models
source('functions.R')
library(glmnet)
library(MASS)
library(GM)
library(DSfdr)
library(knockoff)

set.seed(1482)

##### Input settings
signal = c(9, 9.5, 10, 10.5, 11)
n = 600; p = 1000; s = 20; FDR = 0.1; rho = 0.1
A = signal[1]; beta = rep(0,p)
signal_index = sample((1:p), s, replace=F)
beta[signal_index]= A*sqrt(log(p)/n) 


## Design matrix correlations
sig = toeplitz(seq(rho, 0, length.out = p/10))
V = bdiag(rep(list(sig), 10) ) + diag(rep(1-rho,p))
V = as.matrix(0.1*V)



trials = 100
num_core = 5
cl <- parallel::makeCluster(num_core, setup_strategy = "sequential")
registerDoParallel(cl)
ptime <- system.time({
  r1 <- foreach(icount(trials), .combine=cbind) %dopar% {
    library(glmnet)
    library(MASS)
    library(GM)
    library(DSfdr)
    library(knockoff)
    
    # Generate data 
    X = mvrnorm(n, mu = rep(0,p), Sigma = V)
    #X = scale(X)
    prob = exp(X %*% beta)/(1+exp(X %*% beta))
    y = rep(1,n)
    while(sum(y)/n<0.05 | sum(y)/n>0.95 ){
      for(gen.y in 1:n){
        y[gen.y]=rbinom(1,1,prob[gen.y])
      }
    }
    
    # Implement each method 
    PGM = AdaFDR(x = X, y = y, family = "binomial", kappa = seq(1,10, 1),
                sigma = 1, fdr = FDR, off = 1, reduce = 0)
    
    ds = DS_glm(x=X, y=y, signal_index = signal_index, q = FDR)
    mds =MDS_glm(x=X, y=y, signal_index = signal_index, num_split = 30, q = FDR)
    
    
    k_stat = function(X, Xk, y) stat.lasso_lambdadiff_bin(X, Xk, y, nlambda=200)
    knockoffs = function(X) create.gaussian(X, mu = rep(0,p), Sigma = V)
    
    MX_Knock = knockoff.filter(X, y, knockoffs = knockoffs, statistic=k_stat, fdr = FDR, offset = 1)

    fp = fdp_power(selected_index = PGM$selected, signal_index = signal_index)
    fp_MX = fdp_power(selected_index = MX_Knock$selected, signal_index = signal_index)
    
    c(fp$fdp, fp$power,
      fp_MX$fdp, fp_MX$power,
      ds$DS_fdp, ds$DS_power,
      mds$MDS_fdp, mds$MDS_power
    )
  }
})[3]
stopCluster(cl)
summary_result = round(rowMeans(r1), 3)
data_save1 <- data.frame(Ada_fdp  = summary_result[1],  
                         Ada_power  = summary_result[2],
                         MX_fdp = summary_result[3],
                         MX_power = summary_result[4],
                         DS_fdp = summary_result[5],
                         DS_power = summary_result[6],
                         MDS_fdp = summary_result[7],
                         MDS_power = summary_result[8]
)
data_save1


## Simulations under different signals
set.seed(1482)

A = signal[2]
beta[signal_index]= A*sqrt(log(p)/n) 

cl <- parallel::makeCluster(num_core, setup_strategy = "sequential")
registerDoParallel(cl)
ptime <- system.time({
  r2 <- foreach(icount(trials), .combine=cbind) %dopar% {
    library(glmnet)
    library(MASS)
    library(GM)
    library(DSfdr)
    library(knockoff)
    
    X = mvrnorm(n, mu = rep(0,p), Sigma = V)
    
    prob = exp(X %*% beta)/(1+exp(X %*% beta))
    y = rep(1,n)
    while(sum(y)/n<0.05 | sum(y)/n>0.95 ){
      for(gen.y in 1:n){
        y[gen.y]=rbinom(1,1,prob[gen.y])
      }
    }
    
    PGM = AdaFDR(x = X, y = y, family = "binomial", kappa = seq(1,10, 1),
                sigma = 1, fdr = FDR, off = 1, reduce = 0)
    
  
    ds = DS_glm(x=X, y=y, signal_index = signal_index, q = FDR)
    mds =MDS_glm(x=X, y=y, signal_index = signal_index, num_split = 30, q = FDR)
    
    k_stat = function(X, Xk, y) stat.lasso_lambdadiff_bin(X, Xk, y, nlambda=200)
    knockoffs = function(X) create.gaussian(X, mu = rep(0,p), Sigma = V)
    
    MX_Knock = knockoff.filter(X, y, knockoffs = knockoffs, statistic=k_stat, fdr = FDR, offset = 1)
    
    fp = fdp_power(selected_index = PGM$selected, signal_index = signal_index)
    fp_MX = fdp_power(selected_index = MX_Knock$selected, signal_index = signal_index)
    
    c(fp$fdp, fp$power,
      fp_MX$fdp, fp_MX$power,
      ds$DS_fdp, ds$DS_power,
      mds$MDS_fdp, mds$MDS_power
    )
  }
})[3]
stopCluster(cl)
summary_result = round(rowMeans(r2), 3)
data_save2 <- data.frame(Ada_fdp  = summary_result[1],  
                         Ada_power  = summary_result[2],
                         MX_fdp = summary_result[3],
                         MX_power = summary_result[4],
                         DS_fdp = summary_result[5],
                         DS_power = summary_result[6],
                         MDS_fdp = summary_result[7],
                         MDS_power = summary_result[8]
)
data_save2





set.seed(1482)

A = signal[3]
beta[signal_index]= A*sqrt(log(p)/n) 

cl <- parallel::makeCluster(num_core, setup_strategy = "sequential")
registerDoParallel(cl)
ptime <- system.time({
  r3 <- foreach(icount(trials), .combine=cbind) %dopar% {
    library(glmnet)
    library(MASS)
    library(GM)
    library(DSfdr)
    library(knockoff)
    
    X = mvrnorm(n, mu = rep(0,p), Sigma = V)
    prob = exp(X %*% beta)/(1+exp(X %*% beta))
    y = rep(1,n)
    while(sum(y)/n<0.05 | sum(y)/n>0.95 ){
      for(gen.y in 1:n){
        y[gen.y]=rbinom(1,1,prob[gen.y])
      }
    }
    
    PGM = AdaFDR(x = X, y = y, family = "binomial", kappa = seq(1,10, 1),
                sigma = 1, fdr = FDR, off = 1, reduce = 0)
    
    ds = DS_glm(x=X, y=y, signal_index = signal_index, q = FDR)
    mds =MDS_glm(x=X, y=y, signal_index = signal_index, num_split = 30, q = FDR)
    
    
    k_stat = function(X, Xk, y) stat.lasso_lambdadiff_bin(X, Xk, y, nlambda=200)
    knockoffs = function(X) create.gaussian(X, mu = rep(0,p), Sigma = V)
    
    MX_Knock = knockoff.filter(X, y, knockoffs = knockoffs, statistic=k_stat, fdr = FDR, offset = 1)
    
    fp = fdp_power(selected_index = PGM$selected, signal_index = signal_index)
    
    fp_MX = fdp_power(selected_index = MX_Knock$selected, signal_index = signal_index)
    
    c(fp$fdp, fp$power,
      fp_MX$fdp, fp_MX$power,
      ds$DS_fdp, ds$DS_power,
      mds$MDS_fdp, mds$MDS_power
    )
  }
})[3]
stopCluster(cl)
summary_result = round(rowMeans(r3), 3)
data_save3 <- data.frame(Ada_fdp  = summary_result[1],  
                         Ada_power  = summary_result[2],
                         MX_fdp = summary_result[3],
                         MX_power = summary_result[4],
                         DS_fdp = summary_result[5],
                         DS_power = summary_result[6],
                         MDS_fdp = summary_result[7],
                         MDS_power = summary_result[8]
)
data_save3



set.seed(482)

A = signal[4]
beta[signal_index]= A*sqrt(log(p)/n) 

cl <- parallel::makeCluster(num_core, setup_strategy = "sequential")
registerDoParallel(cl)
ptime <- system.time({
  r4 <- foreach(icount(trials), .combine=cbind) %dopar% {
    library(glmnet)
    library(MASS)
    library(GM)
    library(DSfdr)
    library(knockoff)
    
    X = mvrnorm(n, mu = rep(0,p), Sigma = V)
    prob = exp(X %*% beta)/(1+exp(X %*% beta))
    y = rep(1,n)
    while(sum(y)/n<0.05 | sum(y)/n>0.95 ){
      for(gen.y in 1:n){
        y[gen.y]=rbinom(1,1,prob[gen.y])
      }
    }
    
    PGM = AdaFDR(x = X, y = y, family = "binomial", kappa = seq(1,10, 1),
                sigma = 1, fdr = FDR, off = 1, reduce = 0)

    ds = DS_glm(x=X, y=y, signal_index = signal_index, q = FDR)
    mds =MDS_glm(x=X, y=y, signal_index = signal_index, num_split = 30, q = FDR)
    
    k_stat = function(X, Xk, y) stat.lasso_lambdadiff_bin(X, Xk, y, nlambda=200)
    knockoffs = function(X) create.gaussian(X, mu = rep(0,p), Sigma = V)
   
    MX_Knock = knockoff.filter(X, y, knockoffs = knockoffs, statistic=k_stat, fdr = FDR, offset = 1)
    
    fp = fdp_power(selected_index = PGM$selected, signal_index = signal_index)
    fp_MX = fdp_power(selected_index = MX_Knock$selected, signal_index = signal_index)
    
    c(fp$fdp, fp$power,
      fp_MX$fdp, fp_MX$power,
      ds$DS_fdp, ds$DS_power,
      mds$MDS_fdp, mds$MDS_power
    )
  }
})[3]
stopCluster(cl)
summary_result = round(rowMeans(r4), 3)
data_save4 <- data.frame(Ada_fdp  = summary_result[1],  
                         Ada_power  = summary_result[2],
                         MX_fdp = summary_result[3],
                         MX_power = summary_result[4],
                         DS_fdp = summary_result[5],
                         DS_power = summary_result[6],
                         MDS_fdp = summary_result[7],
                         MDS_power = summary_result[8]
)
data_save4




set.seed(1482)

A = signal[5]
beta[signal_index]= A*sqrt(log(p)/n) 

cl <- parallel::makeCluster(num_core, setup_strategy = "sequential")
registerDoParallel(cl)
ptime <- system.time({
  r5 <- foreach(icount(trials), .combine=cbind) %dopar% {
    library(glmnet)
    library(MASS)
    library(GM)
    library(DSfdr)
    library(knockoff)
    
    X = mvrnorm(n, mu = rep(0,p), Sigma = V)
    prob = exp(X %*% beta)/(1+exp(X %*% beta))
    y = rep(1,n)
    while(sum(y)/n<0.05 | sum(y)/n>0.95 ){
      for(gen.y in 1:n){
        y[gen.y]=rbinom(1,1,prob[gen.y])
      }
    }
    
    PGM = AdaFDR(x = X, y = y, family = "binomial", kappa = seq(1,10, 1),
                sigma = 1, fdr = FDR, off = 1, reduce = 0)
  
    ds = DS_glm(x=X, y=y, signal_index = signal_index, q = FDR)
    mds =MDS_glm(x=X, y=y, signal_index = signal_index, num_split = 30, q = FDR)
    
    k_stat = function(X, Xk, y) stat.lasso_lambdadiff_bin(X, Xk, y, nlambda=200)
    knockoffs = function(X) create.gaussian(X, mu = rep(0,p), Sigma = V)
    
    MX_Knock = knockoff.filter(X, y, knockoffs = knockoffs, statistic=k_stat, fdr = FDR, offset = 1)
    
    fp = fdp_power(selected_index = PGM$selected, signal_index = signal_index)
    
    fp_MX = fdp_power(selected_index = MX_Knock$selected, signal_index = signal_index)
    
    c(fp$fdp, fp$power,
      fp_MX$fdp, fp_MX$power,
      ds$DS_fdp, ds$DS_power,
      mds$MDS_fdp, mds$MDS_power
    )
  }
})[3]
stopCluster(cl)
summary_result = round(rowMeans(r5), 3)
data_save5 <- data.frame(Ada_fdp  = summary_result[1],  
                         Ada_power  = summary_result[2],
                         MX_fdp = summary_result[3],
                         MX_power = summary_result[4],
                         DS_fdp = summary_result[5],
                         DS_power = summary_result[6],
                         MDS_fdp = summary_result[7],
                         MDS_power = summary_result[8]
)
data_save5


mydata = rbind(data_save1, data_save2, 
              data_save3, data_save4, data_save5)
mydata



write.csv(mydata, file = 'logistic_model_final.csv')



mydata = read.csv('logistic_model_final.csv')[-1]


method = c(rep('BGM', 5), rep('Knockoff', 5), 
            rep('DS', 5), rep('MDS', 5))

Signal = rep(signal, 2)
data.frame(Signal, method)
fdrs = c(mydata[,1], mydata[,3], mydata[,5], mydata[,7])
powers = c(mydata[,2], mydata[,4], mydata[,6], mydata[,8])

### calculate results
my_fdr = data.frame(Signal, method, fdrs)
my_power = data.frame(Signal, method, powers)


### plots
library(ggplot2)
library(gridExtra)

p_fdr = ggplot(my_fdr, aes(x = Signal, y = fdrs, col = method)) + 
				geom_point(aes(x = Signal, y = fdrs, shape = method), data = my_fdr, alpha = 1.5, size = 3.5) +
				ylab('FDRs') + 
        geom_line(aes(x = Signal, y = fdrs, linetype = method), data = my_fdr, alpha = 1.75, size = 1.25) + 
				xlab('Signal amplitude') +
        theme_bw() + theme(legend.position=c(0.12,0.8)) + 
        coord_cartesian(ylim = c(0.0, 0.4)) +
        scale_y_continuous(breaks = seq(0,1, 0.1)) +
        geom_hline(aes(yintercept = FDR), colour="#990000", linetype="dashed") +
        theme(legend.title = element_text(size=20),legend.text=element_text(size=18))+
        theme(axis.title.x = element_text(size=25),axis.title.y = element_text(size=25),axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))


p_power = ggplot(my_power, aes(x = Signal, y = powers, col = method)) + 
				geom_point(aes(x = Signal, y = powers, shape = method), data = my_power, alpha = 1.5, size = 5) +
				ylab('Powers') + 
        geom_line(aes(x = Signal, y = powers, linetype = method), data = my_power, alpha = 1.75, size = 1.25) + 
				xlab('Signal amplitude') +
        theme_bw() + theme(legend.position = "none") + 
        coord_cartesian(ylim = c(0.0, 1))  +
        theme(axis.title.x = element_text(size=25),axis.title.y = element_text(size=25),axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
grid.arrange(p_fdr, p_power, ncol=2)

ggsave("logistic_model_final.pdf", arrangeGrob(p_fdr, p_power, ncol=2), width = 18, height = 8)

