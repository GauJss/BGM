
### simulation results with linear models
setwd("/Users/wuyujia/Desktop/KnockOff/Code0227/simulations")
source('functions.R')
library(glmnet)
library(MASS)
library(GM)
library(DSfdr)
library(knockoff)


##### Input settings
set.seed(871)
n = 400; p = 1000; s = 50; FDR = 0.1; rho = 0.2
A = 3.0; beta = rep(0,p)
signal_index = sample((1:p), s, replace=F)
beta[signal_index]= round(runif(s, A*sqrt(log(p)/n), A*sqrt(log(p)/n)+0.2), 3)

## Design matrix correlations
sig = toeplitz(seq(rho, 0, length.out = p/10))
V = bdiag(rep(list(sig), 10) ) + diag(rep(1-rho,p))
V = as.matrix(V)

trials = 50
cl <- parallel::makeCluster(5, setup_strategy = "sequential")
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
    y = X%*%beta + rnorm(n, 0, 1)
    
    # Implement each method  
    PGM = AdaFDR(x = X, y = y, family = "gaussian",  
                sigma = 1, fdr = FDR, off = 1, reduce = 0)
    
    GM = gm(x = X, y = y, q = FDR)
    
    DS=DS_lm(x=X, y=y, q = FDR)
    MDS=MDS_lm(x=X, y=y, num_split=30, q = FDR)
    
    k_stat = function(X, Xk, y) stat.glmnet_coefdiff(X, Xk, y, nfolds=5)
    MX_Knock = knockoff.filter(X, y, knockoffs = create.second_order, statistic=k_stat, fdr = FDR, offset = 1)

    fp = fdp_power(selected_index = PGM$selected, signal_index = signal_index)
    fp_GM = fdp_power(selected_index = GM$gm_selected, signal_index = signal_index)
    fp_DS = fdp_power(selected_index = DS$select_index, signal_index = signal_index)
    fp_MDS = fdp_power(selected_index = MDS$select_index, signal_index = signal_index)
    fp_MX = fdp_power(selected_index = MX_Knock$selected, signal_index = signal_index)
    
    c(fp$fdp, fp$power,
      fp_GM$fdp, fp_GM$power,
      fp_DS$fdp, fp_DS$power,
      fp_MDS$fdp, fp_MDS$power,
      fp_MX$fdp, fp_MX$power
    )
  }
})[3]
stopCluster(cl)
summary_result = round(rowMeans(r1), 3)
data_save1 <- data.frame(Ada_fdp  = summary_result[1],  
                         Ada_power  = summary_result[2],
                         GM_fdp  = summary_result[3],  
                         GM_power  = summary_result[4],
                         DS_fdp  = summary_result[5],  
                         DS_power  = summary_result[6],
                         MDS_fdp  = summary_result[7],  
                         MDS_power  = summary_result[8],
                         MX_fdp = summary_result[9],
                         MX_power = summary_result[10]
)
data_save1


## Simulations under different correlations
set.seed(871)
rho = 0.3
sig = toeplitz(seq(rho, 0, length.out = p/10))
V = bdiag(rep(list(sig), 10) ) + diag(rep(1-rho,p))
V = as.matrix(V)

cl <- parallel::makeCluster(5, setup_strategy = "sequential")
registerDoParallel(cl)
ptime <- system.time({
  r2 <- foreach(icount(trials), .combine=cbind) %dopar% {
    library(glmnet)
    library(MASS)
    library(GM)
    library(DSfdr)
    library(knockoff)
    
    X = mvrnorm(n, mu = rep(0,p), Sigma = V)
    y = X%*%beta + rnorm(n, 0, 1)
    
    PGM = AdaFDR(x = X, y = y, family = "gaussian", 
                sigma = 1, fdr = FDR, off = 1, reduce = 0)
    
    GM = gm(x = X, y = y, q = FDR)
    
    DS=DS_lm(x=X, y=y, q = FDR)
    MDS=MDS_lm(x=X, y=y, num_split=30, q = FDR)
    
    
    k_stat = function(X, Xk, y) stat.glmnet_coefdiff(X, Xk, y, nfolds=5)
    MX_Knock = knockoff.filter(X, y, knockoffs = create.second_order, statistic=k_stat, fdr = FDR, offset = 1)
    
    fp = fdp_power(selected_index = PGM$selected, signal_index = signal_index)
    fp_GM = fdp_power(selected_index = GM$gm_selected, signal_index = signal_index)
    fp_DS = fdp_power(selected_index = DS$select_index, signal_index = signal_index)
    fp_MDS = fdp_power(selected_index = MDS$select_index, signal_index = signal_index)
    fp_MX = fdp_power(selected_index = MX_Knock$selected, signal_index = signal_index)
    
    
    c(fp$fdp, fp$power,
      fp_GM$fdp, fp_GM$power,
      fp_DS$fdp, fp_DS$power,
      fp_MDS$fdp, fp_MDS$power,
      fp_MX$fdp, fp_MX$power
    )
  }
})[3]
stopCluster(cl)
summary_result = round(rowMeans(r2), 3)
data_save2 <- data.frame(Ada_fdp  = summary_result[1],  
                         Ada_power  = summary_result[2],
                         GM_fdp  = summary_result[3],  
                         GM_power  = summary_result[4],
                         DS_fdp  = summary_result[5],  
                         DS_power  = summary_result[6],
                         MDS_fdp  = summary_result[7],  
                         MDS_power  = summary_result[8],
                         MX_fdp = summary_result[9],
                         MX_power = summary_result[10]
)
data_save2




set.seed(871)
rho = 0.4
sig = toeplitz(seq(rho, 0, length.out = p/10))
V = bdiag(rep(list(sig), 10) ) + diag(rep(1-rho,p))
V = as.matrix(V)

cl <- parallel::makeCluster(5, setup_strategy = "sequential")
registerDoParallel(cl)
ptime <- system.time({
  r3 <- foreach(icount(trials), .combine=cbind) %dopar% {
    library(glmnet)
    library(MASS)
    library(GM)
    library(DSfdr)
    library(knockoff)
    
    X = mvrnorm(n, mu = rep(0,p), Sigma = V)
    y = X%*%beta + rnorm(n, 0, 1)
    
    PGM = AdaFDR(x = X, y = y, family = "gaussian",  
                sigma = 1, fdr = FDR, off = 1, reduce = 0)
    
    GM = gm(x = X, y = y, q = FDR)
    
    DS=DS_lm(x=X, y=y, q = FDR)
    MDS=MDS_lm(x=X, y=y, num_split=30, q = FDR)
    
    
    k_stat = function(X, Xk, y) stat.glmnet_coefdiff(X, Xk, y, nfolds=5)
    MX_Knock = knockoff.filter(X, y, knockoffs = create.second_order, statistic=k_stat, fdr = FDR, offset = 1)
    
    fp = fdp_power(selected_index = PGM$selected, signal_index = signal_index)
    fp_GM = fdp_power(selected_index = GM$gm_selected, signal_index = signal_index)
    fp_DS = fdp_power(selected_index = DS$select_index, signal_index = signal_index)
    fp_MDS = fdp_power(selected_index = MDS$select_index, signal_index = signal_index)
    fp_MX = fdp_power(selected_index = MX_Knock$selected, signal_index = signal_index)
    
    c(fp$fdp, fp$power,
      fp_GM$fdp, fp_GM$power,
      fp_DS$fdp, fp_DS$power,
      fp_MDS$fdp, fp_MDS$power,
      fp_MX$fdp, fp_MX$power
    )
  }
})[3]
stopCluster(cl)
summary_result = round(rowMeans(r3), 3)
data_save3 <- data.frame(Ada_fdp  = summary_result[1],  
                         Ada_power  = summary_result[2],
                         GM_fdp  = summary_result[3],  
                         GM_power  = summary_result[4],
                         DS_fdp  = summary_result[5],  
                         DS_power  = summary_result[6],
                         MDS_fdp  = summary_result[7],  
                         MDS_power  = summary_result[8],
                         MX_fdp = summary_result[9],
                         MX_power = summary_result[10]
)
data_save3



set.seed(871)
rho = 0.5
sig = toeplitz(seq(rho, 0, length.out = p/10))
V = bdiag(rep(list(sig), 10) ) + diag(rep(1-rho,p))
V = as.matrix(V)

cl <- parallel::makeCluster(5, setup_strategy = "sequential")
registerDoParallel(cl)
ptime <- system.time({
  r4 <- foreach(icount(trials), .combine=cbind) %dopar% {
    library(glmnet)
    library(MASS)
    library(GM)
    library(DSfdr)
    library(knockoff)
    
    X = mvrnorm(n, mu = rep(0,p), Sigma = V)
    y = X%*%beta + rnorm(n, 0, 1)
    
    PGM = AdaFDR(x = X, y = y, family = "gaussian",  
                sigma = 1, fdr = FDR, off = 1, reduce = 0)
    
    GM = gm(x = X, y = y, q = FDR)
    
    DS=DS_lm(x=X, y=y, q = FDR)
    MDS=MDS_lm(x=X, y=y, num_split=30, q = FDR)
    
    
    k_stat = function(X, Xk, y) stat.glmnet_coefdiff(X, Xk, y, nfolds=5)
    MX_Knock = knockoff.filter(X, y, knockoffs = create.second_order, statistic=k_stat, fdr = FDR, offset = 1)
    
    fp = fdp_power(selected_index = PGM$selected, signal_index = signal_index)
    fp_GM = fdp_power(selected_index = GM$gm_selected, signal_index = signal_index)
    fp_DS = fdp_power(selected_index = DS$select_index, signal_index = signal_index)
    fp_MDS = fdp_power(selected_index = MDS$select_index, signal_index = signal_index)
    fp_MX = fdp_power(selected_index = MX_Knock$selected, signal_index = signal_index)
    
    c(fp$fdp, fp$power,
      fp_GM$fdp, fp_GM$power,
      fp_DS$fdp, fp_DS$power,
      fp_MDS$fdp, fp_MDS$power,
      fp_MX$fdp, fp_MX$power
    )
  }
})[3]
stopCluster(cl)
summary_result = round(rowMeans(r4), 3)
data_save4 <- data.frame(Ada_fdp  = summary_result[1],  
                         Ada_power  = summary_result[2],
                         GM_fdp  = summary_result[3],  
                         GM_power  = summary_result[4],
                         DS_fdp  = summary_result[5],  
                         DS_power  = summary_result[6],
                         MDS_fdp  = summary_result[7],  
                         MDS_power  = summary_result[8],
                         MX_fdp = summary_result[9],
                         MX_power = summary_result[10]
)
data_save4



set.seed(871)
rho = 0.6
sig = toeplitz(seq(rho, 0, length.out = p/10))
V = bdiag(rep(list(sig), 10) ) + diag(rep(1-rho,p))
V = as.matrix(V)

cl <- parallel::makeCluster(5, setup_strategy = "sequential")
registerDoParallel(cl)
ptime <- system.time({
  r5 <- foreach(icount(trials), .combine=cbind) %dopar% {
    library(glmnet)
    library(MASS)
    library(GM)
    library(DSfdr)
    library(knockoff)
    
    X = mvrnorm(n, mu = rep(0,p), Sigma = V)
    y = X%*%beta + rnorm(n, 0, 1)
    
    PGM = AdaFDR(x = X, y = y, family = "gaussian",  
                sigma = 1, fdr = FDR, off = 1, reduce = 0)
    
    GM = gm(x = X, y = y, q = FDR)
    
    DS=DS_lm(x=X, y=y, q = FDR)
    MDS=MDS_lm(x=X, y=y, num_split=30, q = FDR)
    
    
    k_stat = function(X, Xk, y) stat.glmnet_coefdiff(X, Xk, y, nfolds=5)
    MX_Knock = knockoff.filter(X, y, knockoffs = create.second_order, statistic=k_stat, fdr = FDR, offset = 1)
    
    fp = fdp_power(selected_index = PGM$selected, signal_index = signal_index)
    fp_GM = fdp_power(selected_index = GM$gm_selected, signal_index = signal_index)
    fp_DS = fdp_power(selected_index = DS$select_index, signal_index = signal_index)
    fp_MDS = fdp_power(selected_index = MDS$select_index, signal_index = signal_index)
    fp_MX = fdp_power(selected_index = MX_Knock$selected, signal_index = signal_index)
    
    c(fp$fdp, fp$power,
      fp_GM$fdp, fp_GM$power,
      fp_DS$fdp, fp_DS$power,
      fp_MDS$fdp, fp_MDS$power,
      fp_MX$fdp, fp_MX$power
    )
  }
})[3]
stopCluster(cl)
summary_result = round(rowMeans(r5), 3)
data_save5 <- data.frame(Ada_fdp  = summary_result[1],  
                         Ada_power  = summary_result[2],
                         GM_fdp  = summary_result[3],  
                         GM_power  = summary_result[4],
                         DS_fdp  = summary_result[5],  
                         DS_power  = summary_result[6],
                         MDS_fdp  = summary_result[7],  
                         MDS_power  = summary_result[8],
                         MX_fdp = summary_result[9],
                         MX_power = summary_result[10]
)
data_save5


mydata = rbind(data_save1, data_save2, 
              data_save3, data_save4, data_save5)
mydata


write.csv(mydata, file = 'HD_linear_final.csv')




mydata = read.csv('HD_linear_final.csv')[-1]


method = c(rep('BGM', 5), rep('GM', 5), rep('DS', 5), 
           rep('MDS', 5), rep('Knockoff', 5))
corr = c(0.2,0.3, 0.4, 0.5,0.6)
correlation = rep(corr, 5)
fdrs = c(mydata[,1], mydata[,3], mydata[,5], mydata[,7], mydata[,9])
powers = c(mydata[,2], mydata[,4], mydata[,6], mydata[,8], mydata[,10])

### calculate results
my_fdr = data.frame(correlation, method, fdrs)
my_power = data.frame(correlation, method, powers)


### plots
library(ggplot2)
library(gridExtra)

p_fdr = ggplot(my_fdr, aes(x = correlation, y = fdrs, col = method)) + 
				geom_point(aes(x = correlation, y = fdrs, shape = method), data = my_fdr, alpha = 1.5, size = 3.5) +
				ylab('FDRs') + 
        geom_line(aes(x = correlation, y = fdrs, linetype = method), data = my_fdr, alpha = 1.75, size = 1.25) + 
				xlab('Correlation factor') +
        theme_bw() + theme(legend.position=c(0.12,0.8)) + 
        coord_cartesian(ylim = c(0.0, 0.4)) +
        scale_y_continuous(breaks = seq(0,1, 0.1)) +
        geom_hline(aes(yintercept = FDR), colour="#990000", linetype="dashed") +
        theme(legend.title = element_text(size=20),legend.text=element_text(size=18))+
        theme(axis.title.x = element_text(size=25),axis.title.y = element_text(size=25),axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))


p_power = ggplot(my_power, aes(x = correlation, y = powers, col = method)) + 
				geom_point(aes(x = correlation, y = powers, shape = method), data = my_power, alpha = 1.5, size = 5) +
				ylab('Powers') + 
        geom_line(aes(x = correlation, y = powers, linetype = method), data = my_power, alpha = 1.75, size = 1.25) + 
				xlab('Correlation factor') +
        theme_bw() + theme(legend.position = "none") + 
        coord_cartesian(ylim = c(0.3, 1))  + scale_y_continuous(breaks = seq(0,1, 0.1)) +
        theme(axis.title.x = element_text(size=25),axis.title.y = element_text(size=25),axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
grid.arrange(p_fdr, p_power, ncol=2)

ggsave("HD_linear_final.pdf", arrangeGrob(p_fdr, p_power, ncol=2), width = 18, height = 8)

