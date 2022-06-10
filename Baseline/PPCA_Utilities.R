split <- function(masked_data, len_masked){
  
  # Split Train and Validation Data in terms of the mask length 
  
  train <- masked_data[1:(nrow(masked_data) - len_masked), ]
  
  valid <- masked_data[(nrow(masked_data) - len_masked + 1):nrow(masked_data), ]
  
  list(train=train, valid=valid)
} 

fit_ppca <- function(train, Pcs=10, method="ppca"){
  
  # PCAs
  ppca_ <-  pca(train, nPcs=Pcs, method="ppca")
  
  # Latent Varible aka Sigma in MVNorm
  pc_scores <- scores(ppca_)
  
  # Contribution of each observation to each the Latent variable
  pc_loadings <- loadings(ppca_)
  
  # \mu vector. Mass Center
  pc_mu <- center(ppca_)
  
  # Estimated Noise
  epsilon <- (fitted(ppca_)-train)
  
  train_reconstructed <- fitted(ppca_)
  
  # Estimated Noise Variance
  var_eps <- var(c(epsilon))
  
  fittedval <- fitted(ppca_)
  
  return(list(pc_scores = pc_scores,pc_loadings = pc_loadings, pc_mu = pc_mu, epsilon = epsilon, var_eps = var_eps, fitted =fittedval))
  
}

chol_inv <- function(A){
  chol2inv(chol(A))
}


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
  
  # chol2inv(chol(M)) where M is cov matrix 
  
  # Compute x_mis|x_obs
  # Tolerance for the solve
  sample[mis_ind] <- mu_mis + Sigma_mis %*% chol_inv(t(Sigma_obs) %*% Sigma_obs + var_eps*diag(ncol(Sigma_obs))) %*%  t(Sigma_obs) %*%  (sample_obs - mu_obs) 
  
  #Sigma_cond <- var_eps * Sigma_mis %*% solve(t(Sigma_obs) %*% Sigma_obs + var_eps*diag(ncol(Sigma_obs))) %*% t(Sigma_mis) + var_eps*diag(nrow(Sigma_mis))
  sample
  
}


pc.predict <- function(valid_set, mu, Sigma, var_eps){
  
  # Run pc.predict_1_sample for all functions
  
  val_pred <- apply(valid_set, 1, pc.predict_1_sample, mu, Sigma, var_eps)
  t(val_pred)
  
}

fit_bpca <- function(train, Pcs=10, method="ppca"){
  
  # PCAs
  bpca_ <-  pca(train, nPcs=Pcs, method="bpca")
  
  # Latent Varible aka Sigma in MVNorm
  pc_scores <- scores(bpca_)
  
  # Contribution of each observation to each the Latent variable
  pc_loadings <- loadings(bpca_)
  
  # \mu vector. Mass Center
  pc_mu <- center(bpca_)
  
  # Estimated Noise
  epsilon <- (fitted(bpca_)-train)
  
  train_reconstructed <- fitted(bpca_)
  
  # Estimated Noise Variance
  var_eps <- var(c(epsilon))
  
  fittedval <- fitted(bpca_)
  
  return(list(pc_scores = pc_scores,pc_loadings = pc_loadings, pc_mu = pc_mu, epsilon = epsilon, var_eps = var_eps, fitted =fittedval))
  
}

metrics <- function(valid_unmask, missing, valid_pred){
  
  corrected_Abs_SE <- sqrt(mean(abs(valid_unmask[missing] - valid_pred[missing])^2) / mean(abs(valid_unmask[missing])^2))
  MAE <- mean(abs(valid_unmask[missing] - valid_pred[missing]))
  RMSE <- sqrt(mean((valid_unmask[missing] - valid_pred[missing])^2))
  MSE <- mean((valid_unmask[missing] - valid_pred[missing])^2)
  
  cat(" NRMSE: ", corrected_Abs_SE)
  cat("\n MAE: ", MAE)
  cat("\n RMSE: ",RMSE) 
  cat("\n MSE:", MSE)
  
  return(c(corrected_Abs_SE, MAE, RMSE, MSE))
  
  }

table_metrics <- function(metrics_vec){
  
table <- matrix(metrics_vec, 1, 4)
rownames(table) <- c("Error Metric")
colnames(table) <- c("NRMSE", "MAE", "RMSE", "MSE")

knitr::kable(table)}


plot_pred_vs_true <- function(valid_unmask, missing, valid_pred, ds_name, col='blue'){
  
  plot(valid_unmask[missing], valid_pred[missing], col= col,
       main=paste("Hold-out Prediction PPCA ", ds_name), ylab= 'Prediction', 
       xlab= 'Background Truth') 
  
  #+points(valid_unmask[missing], valid_pred[missing], col='red' ) 
  
  legend(1,95, legend=c("Predictions", "Post-Processed prediction"),
         col=c("red", "blue"), lty=1:2, cex=0.8)
  
}



GridSearch_ppca <- function(train, itNum){
  
  NRMSE_vec <-  rep(NaN, itNum)

  for (k in 2:itNum){
    #
    print(k)
    ppca_ <- fit_ppca(train, Pcs=k, method="ppca")
    
    pc_scores <- ppca_$pc_scores
    pc_loadings <- ppca_$pc_loadings
    pc_mu <- ppca_$pc_mu
    var_eps <- ppca_$var_eps
    
    valid_pred  <- pc.predict(valid, pc_mu, pc_loadings, var_eps)
    
    NRMSE_vec[k] <- mean((valid_unmask[missing] - valid_pred[missing])^2)
  
  
  # MAE <- mean(abs(valid_unmask[missing] - valid_pred[missing]))
  # RMSE <- sqrt(mean((valid_unmask[missing] - valid_pred[missing])^2))
  # corrected_Abs_SE <- sqrt(mean(abs(valid_unmask[missing] - valid_pred[missing])^2) / mean(abs(valid_unmask[missing])^2))
  
  }
  plot(2:(itNum), NRMSE_vec[2:length(NRMSE_vec)])
  
  cat('Number of Optimal PPCs', (which.min(NRMSE_vec) + 1) )
}



plot_1MPC <- function(data, Real, Nmpc, ds_name){
  
  if (Real){
  
    plot(data[,1], main = paste("plot of 1st MPC in 4 subbands Evolution Over Time (Real Part) ", ds_name), xlab = "TTI", ylab="Amplitude", ylim = c(-2,2), type = "l" , cex=1.5)
    points(data[,(1+Nmpc*2)], col=2, type = "l" )
    points(data[,(1+Nmpc*4)], col=3, type = "l" )
    points(data[,(1+Nmpc*6)], col=4, type = "l" )
    legend(x = "top", legend=c("1st MPC SB1", "1st MPC SB2", "1st MPC SB3", "1st MPC SB4"),
         col=c( 1,2,3,4), lty = c(1,1), lwd = c(1,1))}
  
  else{
    
    plot(  data[, (1+Nmpc)], main = paste("plot of 1st MPC in 4 subbands Evolution Over Time (Real Part) ", ds_name), xlab = "TTI", ylab="Amplitude", ylim=c(-2,2.5), type="l")
    points(data[,(1+3*Nmpc)], col=2, type="l")
    points(data[,(1+5*Nmpc)], col=3, type="l")
    points(data[,(1+7*Nmpc)], col=4, type="l")
    legend(x = "top", legend=c("1st MPC SB1", "1st MPC SB2", "1st MPC SB3", "1st MPC SB4"),
           col=c( 1,2,3,4), lty = c(1,1), lwd = c(1,1))
    
  }
  
  
}
  


plot_predictions_1TTI <- function(valid_unmask, valid_pred, samplenum, ylim=NULL, ds_name){
  
  
  # Missing indexes, indexes pred False True
  missingidx <- is.na(valid[samplenum,])
  # Missing indexes, just indexes
  idxnan <- which(missingidx)
  
  # ypred values
  ypred <- valid_pred[samplenum, missingidx]
  
  plot(valid_unmask[samplenum,],  col=1,
       main = paste("Prediction of TTI n =", samplenum, ' ',ds_name),
       ylab = "Amplitude", xlab = "Feature")
  
  points(idxnan, ypred, col=2)
  
  legend(x = "bottomright", legend=c("Predictions", "Background Truth"),
         col=c("red", "black"), lty=1, lwd = 1)
  
}




