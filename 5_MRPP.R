
library(vegan)
library(corrplot)


setwd('/home/javierlopatin/Documents/temp/WetlandPhenology/')

df = read.csv('data/All_variables_random.csv')
colnames(df)
df$PFT = as.factor(df$PFT)
df$cluster = as.factor(df$cluster)

rmse = read.csv('data/rmse_random.csv')
colnames(rmse)

dem = read.csv('data/dem_random.csv')

# correlation analysis of all variables
all = cbind(df, rmse[, 2:ncol(rmse)], dem$dem)
all <- na.omit(all)

M = cor(all[, c(6:21, 75:ncol(all))])
testRes = cor.mtest(all[, c(6:21, 75:ncol(all))], conf.level = 0.95)
corrplot(M, p.mat = testRes$p, sig.level = 0.05, addrect = 2, type="lower")


# MRPP --------------------------------------------------------------------

# create empty matrices
mat_veg  = matrix(NA, ncol=ncol(all[, c(6:21, 75:ncol(all))]), nrow=2)
colnames(mat_veg) = colnames(all[, c(6:21, 75:ncol(all))])
rownames(mat_veg) = c("A", "P")

mat_clust  = matrix(NA, ncol=ncol(all[, c(6:21, 75:ncol(all))]), nrow=2)
colnames(mat_clust) = colnames(all[, c(6:21, 75:ncol(all))])
rownames(mat_clust) = c("A", "P")

# run MRPP through the variables
MRPP <- function(X, Y, mat, num_cores = 6){
  pb <- txtProgressBar(min = 0, max = ncol(X), style = 3,  width = 50, char = "=") # progress bar
  length = ncol(X)
  for(i in 1:length){
    obj_mrpp = mrpp(dat = X[,i], grouping = Y, parallel = num_cores, distance = "mahalanobis")
    mat[1,i] = obj_mrpp$A
    mat[2,i] = obj_mrpp$Pvalue
    setTxtProgressBar(pb, i)
  }
  close(pb) # close progress bar
  return(mat)
}

# vegetation types
mrpp_veg = MRPP(X=all[, c(6:21, 75:ncol(all))], Y=all$PFT, mat=mat_veg)
mrpp_veg

# phenology types
mrpp_clust = MRPP(all[, c(6:21, 75:ncol(all))], all$cluster, mat_clust)
mrpp_clust

save.image('mrpp.RData')


# plots -------------------------------------------------------------------

# functions to plot the relatve importance using inferno colors
varImp <- function(varImport = mrpp_lsp){
  library("viridis")  
  
  # add coefficients
  colorpal <-  inferno(100)
  
  # variables from ensemble
  imp <- varImport
  
  # matices of varImport
  z1 <- matrix (rep (imp, 100), ncol=100)
  
  # MRFF coefficients
  wl = seq(length(imp))
  #image(wl, seq(0, 100, 1), z1, xlim = c(min(wl)-10, max(wl)+10), xlab=expression(lambda(nm)), col=blueish, ylab="", axes=F, cex.lab = 1.3)
  image(wl, seq(0, 100, 1), z1, xlim = c(min(wl)-10, max(wl)+10), xlab="", col=colorpal  , ylab="", axes=F, cex.lab = 1.3)    
}



svg('figures/mrpp_veg.svg')
varImp(mrpp_veg[1, ])
dev.off()

svg('figures/mrpp_clust.svg')
varImp(mrpp_clust[1, ])
dev.off()


# variable importance with PLS-DA -----------------------------------------

library(caret)

tr_control <- trainControl(method = 'cv', number = 10)

# train model using vegetation clases
PLS_veg <- train(x=all[, c(6:21, 75:ncol(all))], y=all$PFT , method = "pls",  tuneLength=20, 
                 trControl = tr_control, preProcess = c("center", "scale"), 
                 metric = "Kappa", maximize = T, na.action=na.omit)
PLS_veg

plot(PLS_veg$results[, 1], PLS_veg$results[, 3], type='l', ylab='kappa')

# train model using cluster clases
PLS_clust <- train(x=all[, c(6:21, 75:ncol(all))], y=all$cluster , method = "pls",  tuneLength=20, 
                   trControl = tr_control, preProcess = c("center", "scale"), 
                   metric = "Kappa", maximize = T, na.action=na.omit)
PLS_clust

plot(PLS_clust$results[, 1], PLS_clust$results[, 3], type='l', ylab='kappa')

relPLSimp = function(X, Y, ncomp, iter=100){
  #############################################################
  #
  # Estimate relative importance using PLS-DA by using Kappa 
  # and predicting values in random variables
  #
  # Y = observed classes (vector of Factor class)
  # X = predictors (data frame of numerical variables)
  # iter = number of bootstrap iterations (intener)
  #
  ##############################################################
  
  kappa <- matrix(nrow = iter, ncol = ncol(X))
  colnames(kappa) = colnames(X)
  
  boot <- createResample(Y, times = iter, list = TRUE)
  
  pb <- txtProgressBar(min = 0, max = iter, style = 3, width = 50, char = "=") 
  
  for (i in 1:iter){ # bootstraps
    
    train_X <- X[boot[[i]],]
    train_Y <- Y[boot[[i]]]  
    val_X   <- X[-boot[[i]],]
    val_Y   <- Y[-boot[[i]]]
    
    PLS  <- caret::plsda(x = train_X, y = train_Y, ncomp = ncomp, probMethod = 'softmax')
    pred <- predict(PLS, val_X)
    
    mat_pls <- confusionMatrix(pred, val_Y)
    
    for (k in 1:(ncol(X))){ # loop through variables
      validar_rand  <- val_X
      kappp = c()
      for (j in 1:10){ # 10 random replaces per variable
        validar_rand[,k]  <- runif(nrow(validar_rand), min=min(validar_rand[,k]), max=max(validar_rand[,k]))
        pred_rand <- predict(PLS, validar_rand)
        mat_pls_rand <- confusionMatrix(pred_rand, val_Y)
        kappp[j] <- mat_pls_rand$overall[2]
      }
      # delta kappa
      kappa[i,k] <- mat_pls$overall[2] - median(kappp)
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(kappa)
}

# process PLSDA variable importance
imp_veg = relPLSimp(X=all[, c(6:21, 75:ncol(all))], Y=all$PFT, ncomp=15, iter = 100)
colMeans(imp_veg)
barplot(colMeans(imp_veg))

svg('pls_veg.svg')
varImp(colMeans(imp_veg))
dev.off()

imp_clust = relPLSimp(X=all[, c(6:21, 75:ncol(all))], Y=all$cluster, ncomp=3, iter = 100)
colMeans(imp_clust)
barplot(colMeans(imp_clust))

svg('pls_cluster.svg')
varImp(colMeans(imp_clust))
dev.off()

save(imp_veg, file='plsImp_veg.RData')
save(imp_clust, file='plsImp_clust.RData')

# load data
load(file='plsImp_veg.RData')
load(file='plsImp_clust.RData')
load('mrpp.RData')

# variable importance
#PLS_veg2 <- update(PLS_veg, param = list(ncomp = 10))
plscf <- as.vector(rowMeans(PLS_veg$finalModel$coefficients)) ## extract coeff.
plscf <- plscf / sd (plscf) ## scale regression coefficients

#svg('Figures/imp_veg_coeff.svg', width = 10, height = 5)
barplot(plscf, las=1)
#dev.off()

# svg('img_veg.svg', width = 12, height = 3)
boxplot(imp_veg, col='gray', outline=FALSE, las=2)
abline(0,0)
# dev.off()

# svg('img_clust.svg', width = 12, height = 3)
boxplot(imp_clust, col='gray', outline=FALSE, las=2)
abline(0,0)
# dev.off()


###########3

rmse = read.csv('/home/javierlopatin/Downloads/points_rmse.csv')
lsp = read.csv('/home/javierlopatin/Downloads/points_LSP.csv')

colnames(rmse) = c("x","rmse_all", "rmse_sos", "rmse_pos", "rmse_eos")
colnames(lsp) = c("x",'SOS', 'POS', 'EOS', 'vSOS', 'vPOS', 'vEOS', 'LOS', 'MSP', 'MAU', 'vMSP', 'vMAU', 'AOS', 'IOS', 'ROG', 'ROS', 'SW')

all = cbind(rmse[,2:5], lsp[, 2:17])
all=na.omit(all)

colnames(all)

# correlation analysis of all variables
M = cor( all )
testRes = cor.mtest( all, conf.level = 0.95)
png("corrplot.png")
corrplot(M, p.mat = testRes$p, sig.level = 0.05, addrect = 2, type="lower")
dev.off()
