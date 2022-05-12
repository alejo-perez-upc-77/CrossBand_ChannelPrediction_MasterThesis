chol_inv <- function(A){
  chol2inv(chol(A))
}

pc.predict_1_sample_varcond <- function(sample, mu, Sigma, var_eps){
  
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
  #sample[mis_ind] <- mu_mis + Sigma_mis %*% solve(t(Sigma_obs) %*% 
  #Sigma_obs + var_eps*diag(ncol(Sigma_obs))) %*%  t(Sigma_obs) %*% 
  #(sample_obs - mu_obs) 
  
  Sigma_cond <- var_eps * Sigma_mis %*% 
    chol_inv(t(Sigma_obs) %*% Sigma_obs
    + var_eps*diag(ncol(Sigma_obs))) %*% t(Sigma_mis) +
    var_eps*diag(nrow(Sigma_mis))
  
  diag_cond <- diag(Sigma_cond)
  dim(diag_cond) <- c( 1, length(diag_cond))
  diag_cond
}

pc.predict_varcond <- function(valid_set, mu, Sigma, var_eps){
  
  # Run pc.predict_1_sample for all functions
  
  var_pred <- apply(valid_set, 1, pc.predict_1_sample_varcond, mu,
                    Sigma, var_eps)
  t(var_pred)
  
}
covdiag_pred <- pc.predict_varcond(valid, pc_mu, pc_loadings, var_eps)

########################################################################
###########  Plot Confidence Interval of the Prediction ################
########################################################################
samplenum = 7

plot_predictions <- function(valid_unmask, valid_pred, samplenum, ylim=NULL){
  
  
  # Missing indexes, indexes pred False True
  missingidx <- is.na(valid[samplenum,])
  # Missing indexes, just indexes
  idxnan <- which(missingidx)
  
  # ypred values
  ypred <- valid_pred[samplenum, missingidx]

  plot(valid_unmask[samplenum,],  col=1,
       main = paste("Prediction of TTI n =", samplenum, " (UE0)"),
       ylab = "Amplitude", xlab = "Feature")
  
  points(idxnan, ypred, col=2)
  
  legend(x = "bottomright", legend=c("Predictions", "Background Truth"),
                col=c("red", "black"), lty=1, lwd = 1)
  
}

missingidx <- is.na(valid[samplenum,])
# Missing indexes, just indexes
idxnan <- which(missingidx)

#ypred <- valid_pred[samplenum, missingidx]



plot_predictions(valid_unmask, valid_pred, samplenum)

#points(idxnan, ypred + 2*sqrt(covdiag_pred[samplenum,]), col=7)


########################################################################
###########  Plot Confidence Interval of the Prediction ################
########################################################################

samplenum = 7

plot_CI <- function(valid, valid_pred, covdiag_pred, samplenum, ylim=NULL){
  

  # Missing indexes, indexes pred False True
  missingidx <- is.na(valid[samplenum,])
  # Missing indexes, just indexes
  idxnan <- which(missingidx)
  
  # ypred values
  ypred <- valid_pred[samplenum, missingidx]
  # sorted y_predvalues
  ypred_sorted <- sort(ypred)
  
  idxord <- order(ypred)
  
  std_pred <- sqrt(covdiag_pred[samplenum,]) * 1.96
  std_pred_sorted <- std_pred[idxord]
  
  #hardcoded
  #xlima <- c(1000, 1500)
  
  plot( ypred_sorted,  col=4, ylim=ylim, main= paste('Sorted Predictions for TTI n = ',samplenum,' + predicted CIs)'))
  
  points(ypred_sorted + std_pred[idxord], col=2)
  
  points(ypred_sorted - std_pred[idxord], col=2)
  
  abline(h=mean(ypred_sorted), lty=2, lwd= 2)
  
  legend(x = "top", legend=c("CI", "Prediction Sorted", "mean(Prediction Sorted)"),
         col=c("red", "blue", "black"), lty = c(1,1,3), lwd = c(1,1,3))
  
  }

plot_CI(valid, valid_pred, covdiag_pred, samplenum, ylim=c(-.2,.2))


#pdf(file='Plot_CI.pdf')
#dev.off()



########################################################################
###########  Plot Reconstruction of one TTI of training ################
########################################################################

samplenum = 30

plot(fitted(ppca_)[samplenum,],  
     main="Plot Reconstruction of one TTI of training")
points(train[samplenum,], col=2)
legend(x = "top", legend=c("Predictions", "Background Truth"),
       col=c("red", "black"), lty=1, lwd = 1)

########################################################################

ord = order(c(valid_unmask[missing]))

plot(1:length(c(valid_unmask[missing])), (c(valid_pred[missing])[ord]), 
     pch=1,col="red", main="Prediction and Background Truth Sorted Values UE0",
     xlab="Feature", ylab="Amplitude")




points(1:length(c(valid_unmask[missing])), sort(c(valid_unmask[missing])))

legend(x = "topleft", legend=c("Prediction", "Background Truth"),
       col=c("red", "black"), lty = c(1,1), lwd = c(1,1))









