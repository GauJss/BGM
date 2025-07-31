rm(list = ls())
source('functions.R')
library(knockoff)
library(glmnet)
library(MASS)
library(GM)
library(DSfdr)
library(dplyr)


#replicate <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#set.seed(replicate)

data <- t(read.table("GSE141834_scRNAseq_rawCounts.txt"))
y <- c(rep(0, 400), rep(1, 2000))


### dimension reduction
freq <- colSums(data != 0)/nrow(data)
remove_index <- which(freq < 0.1)
X <- data[, -remove_index]
std <- apply(X, 2, sd)
gene_index <- order(std, decreasing = T)
p <- 500
n <- nrow(X)
X <- X[, gene_index[1:p]]

### FDR control level 
fdr <- 0.1

# Run the knockoff filter.
  knock.gen = function(x) create.second_order(x,  method='equi')
  result = knockoff.filter(X, y, fdr=fdr, knockoffs=knock.gen, statistic=stat.glmnet_lambdasmax)
  knockoff_selected = names(result$selected)
  
  # Run DS & MDS
  num_split=5
  SDS_selected = colnames(X)[DS(X, y, fdr)$SDS]
  MDS_selected = colnames(X)[MDS(X, y, num_split,fdr)$MDS]
  
  ## Run BGM
  PGM = AdaFDR(x = X, y = y, family = "binomial", sigma = 1, fdr = fdr, off = 1, reduce = 1)
  PGM_selected = colnames(X)[PGM$selected]
 

  
  
selection=list(AdaFDR=PGM_selected, DS = SDS_selected, MDS = MDS_selected, 
       Knockoff = knockoff_selected) #, BHq = bhq_selected)  


##commonly selected features in Dai2023 supported by references
ref_truth=c("SERPINA6","FKBP5","NFKBIA","RPL10","SEMA3C",
            "HSPB1","RBBP7","EIF4EBP1","S100A11","NUPR1","MSX2",
            "LY6E","BLOC1S1","HSPA8","HSPA1A","EEF1A1","BCL6","ATF4",
            "IGFBP4","YWHAQ","DDIT4","IRX2","GATA3-AS1","RBM24",
            "TACSTD2","DSCAM-AS1","C1QBP","SNHG19","PRR15L",
            "RPLP0P6","UHMK1")


##commonly selected features by these methods and Dai2023
comm=list(BGM_insec=intersect(PGM_selected,ref_truth), 
     DS_insec=intersect(SDS_selected,ref_truth), MDS_insec=intersect(MDS_selected,ref_truth), 
     Knockoff_insec=intersect(knockoff_selected,ref_truth)) 


##selected features Dai2023 exclude in these methods
diff1=list(BGM_insec=intersect(PGM_selected,ref_truth), 
          DS_insec=intersect(SDS_selected,ref_truth), 
          MDS_insec=intersect(MDS_selected,ref_truth), 
          Knockoff_insec=intersect(knockoff_selected,ref_truth)) 

#selected features exclude in Dai2023 supported by references
diff2=list(BGM_insec=intersect(ref_truth[1:27],PGM_selected), 
           DS_insec=intersect(ref_truth[1:27],SDS_selected), 
           MDS_insec=intersect(ref_truth[1:27],MDS_selected), 
           Knockoff_insec=intersect(ref_truth[1:27],knockoff_selected)) 






