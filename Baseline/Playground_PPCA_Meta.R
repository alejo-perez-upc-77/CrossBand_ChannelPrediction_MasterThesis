library(pcaMethods)
data(metaboliteData)
mD <- metaboliteData
data(metaboliteDataComplete)
mdComp <- metaboliteDataComplete

mis_idx <- which(is.na(mD))


###########################################################################   
############################## PPCA method ################################ 
###########################################################################  

# Extract PCA, Latent variables and \mu vector from training data

# PPCA method
#prep_md <- prep(mD)

ppca_ <-  ppca(mD, 5)

# Latent Varible aka Sigma in MVNorm
pc_scores <- scores(ppca_)

# Contribution of each observation to each the Latent variable
pc_loadings <- loadings(ppca_)

# \mu vector. Mass Center
pc_mu <- center(ppca_)

product_fit <- pc_scores %*% t(pc_loadings)


compObs <- completeObs(ppca_)


fit <- fitted(ppca_)

#mD[mis_idx] <- imputed[mis_idx]



plot(metaboliteDataComplete[mis_idx], compObs[mis_idx] )


plot(metaboliteDataComplete[mis_idx], fit[mis_idx], col='green' )
plot(metaboliteDataComplete[mis_idx], product_fit[mis_idx] , col='red')

###########################################################################   
######################### PCA + PPCA method ############################### 
########################################################################### 

ppca_ <-  pca(mD,'ppca', 5)

# Latent Varible aka Sigma in MVNorm
pc_scores <- scores(ppca_)

# Contribution of each observation to each the Latent variable
pc_loadings <- loadings(ppca_)

# \mu vector. Mass Center
pc_mu <- center(ppca_)

product_fit <- pc_scores %*% t(pc_loadings)
recData <- prep(product_fit, scl(ppca_), center(ppca_), reverse=TRUE)


compObs <- completeObs(ppca_)


fit <- fitted(ppca_)

#mD[mis_idx] <- imputed[mis_idx]

plot(metaboliteDataComplete[mis_idx], compObs[mis_idx], ylim = c(-2,2.5) )


plot(metaboliteDataComplete[mis_idx], fit[mis_idx], col='green', ylim = c(-2,2.5) )
plot(metaboliteDataComplete[mis_idx], recData[mis_idx] , col='red', 
     ylim = c(-2,2.5))
plot(metaboliteDataComplete[mis_idx], product_fit[mis_idx] , col='purple', 
     ylim = c(-2,2.5))
###########################################################################   
############### PCA + PPCA method more with preprocessing ################# 
###########################################################################

prepres <- prep(mD, scale='none', center=TRUE, simple = FALSE)
res <- ppca(prepres$data, nPcs=5) 
missing <- is.na(mD)
res@missing <- missing
res@center <- prepres$center
res@scale <- prepres$scale

cObs <- mD 
product_fit <- scores(res) %*% t(loadings(res))
recData <- prep(product_fit, scl(res), center(res), reverse=TRUE)

cObs[missing] <- recData[missing]


plot(metaboliteDataComplete[mis_idx], recData[missing], ylim = c(-2,2.5) )


###########################################################################   
######################### Functions Inference ############################# 
###########################################################################  

pc.predict_1_sample <- function(sample, mu, Sigma, var_eps){
  
  # This function predicts the missing data for one sample x_i by means of the 
  # conditional predictive distribution
  
  # Split Observations x_i
  obs_ind <- which(!is.na(sample))
  mis_ind <- which(is.na(sample))
  
  sample_obs <- sample[obs_ind]
  
  # Split Mus
  mu_obs <- mu[obs_ind]
  mu_mis <- mu[mis_ind]
  
  # Split Sigma (W)
  Sigma_obs <- Sigma[obs_ind, ]
  Sigma_mis <- Sigma[mis_ind, ]
  
  
  # Compute x_mis|x_obs
  sample[mis_ind] <- mu_mis + Sigma_mis %*% solve(t(Sigma_obs) %*% Sigma_obs + var_eps*diag(ncol(Sigma_obs))) %*%  t(Sigma_obs) %*%  (sample_obs - mu_obs) 
  
  #Sigma_cond <- var_eps * Sigma_mis %*% solve(t(Sigma_obs) %*% Sigma_obs + var_eps*diag(ncol(Sigma_obs))) %*% t(Sigma_mis) + var_eps*diag(nrow(Sigma_mis))
  sample
  
}

pc.predict <- function(valid_set, mu, Sigma, var_eps){
  
  # Run pc.predict_1_sample for all functions
  
  val_pred <- apply(valid_set, 1, pc.predict_1_sample, mu, Sigma, var_eps)
  t(val_pred)
  
}
###########################################################################   
#################  Inference in validation Metabolite Set #################
########################################################################### 

######################## generate train and test ##########################


dim(mdComp)
dim(mD)
N_train <- 90

train <- mdComp[1:N_train,]
valid <- mD[(N_train+1):nrow(mD),]
valid_true <- mdComp[(N_train+1):nrow(mdComp),]


missing <- is.na(valid)

cat("Percentage of NaN Values in Training Set:", (sum(is.na(valid)))/length(valid)*100, '%')

########################### Extraction of PC #############################

ppca_ <-  pca(train,'ppca', 5)

# Latent Varible aka Sigma in MVNorm
pc_scores <- scores(ppca_)

# Contribution of each observation to each the Latent variable
pc_loadings <- loadings(ppca_)

# \mu vector. Mass Center
pc_mu <- center(ppca_)

fit <- fitted(ppca_)

var_eps <- var(c(fit - train))

rec_data <- pc.predict(valid, pc_mu, pc_loadings, var_eps)

rec_data_postpro <- prep(rec_data, scl(ppca_), center(ppca_), reverse=TRUE)


plot(rec_data[missing], valid_true[missing]) + 
  points(rec_data_postpro[missing], valid_true[missing], col='red')

















