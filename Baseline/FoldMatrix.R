
reIm2complex_sb <-  function(sbi){
  # It takes a matrix of 1 subband but Re and Im part concatenated. It merges 
  # them together into a complex subband of MPCs
  
  # Inputs
  # sbi: matrix [tti x p·2]
  
  # Outputs
  # sbfold: matrix [tti x p]
  

  ncol_out <- ncol(sbi)/2
  nrow_out <- nrow(sbi)
  
  sbfold <- matrix(nrow = nrow_out, ncol = ncol_out)
  
  for (i in 1:ncol_out){
    sbfold[, i] <- complex(real=sbi[,i], imaginary=sbi[,i+ncol_out]) }
  
  sbfold
}



fold_matrix <- function(unfold_mat){
  # It takes an unfold matrix of MPCs and returns a forded one
  
  # Inputs
  # sbi: matrix [tti x p·2·sb]
  
  # Outputs
  # sbfold: matrix [p  x sb x tti]
  
  sb = 4
  sbcolreim <- ncol(unfold_mat)/sb
  nrows_sb <- nrow(unfold_mat)
  
  p <- sbcolreim/2
  
  sb1 <- unfold_mat[, 1:sbcolreim]
  sb2 <- unfold_mat[, (sbcolreim+ 1) : (2*sbcolreim)]
  sb3 <- unfold_mat[, (2*sbcolreim + 1) : (3*sbcolreim)]
  sb4 <- unfold_mat[, (3*sbcolreim + 1) : ncol(unfold_mat)]
  
  sb1_fold <- reIm2complex_sb(sb1) 
  sb2_fold <- reIm2complex_sb(sb2) 
  sb3_fold <- reIm2complex_sb(sb3) 
  sb4_fold <- reIm2complex_sb(sb4) 
  
  fold_mat <-  sb1_fold
  fold_mat <- abind(fold_mat, sb2_fold, along = 3)
  fold_mat <- abind(fold_mat, sb3_fold, along = 3)
  fold_mat <- abind(fold_mat, sb4_fold, along = 3)
  
  aperm(fold_mat, c(2,3,1))
  
}

fold_matrix_halfway <- function(unfold_mat){
  # It takes an unfold matrix of MPCs and returns a forded one
  
  # Inputs
  # sbi: matrix [tti x p·2·sb]
  
  # Outputs
  # sbfold: matrix [p  x sb x tti]
  
  sb = 4
  sbcolreim <- ncol(unfold_mat)/sb
  nrows_sb <- nrow(unfold_mat)
  
  p <- sbcolreim/2
  
  sb1 <- unfold_mat[, 1:sbcolreim]
  sb2 <- unfold_mat[, (sbcolreim+ 1) : (2*sbcolreim)]
  sb3 <- unfold_mat[, (2*sbcolreim + 1) : (3*sbcolreim)]
  sb4 <- unfold_mat[, (3*sbcolreim + 1) : ncol(unfold_mat)]
  
  sb1_fold <- reIm2complex_sb(sb1) 
  sb2_fold <- reIm2complex_sb(sb2) 
  sb3_fold <- reIm2complex_sb(sb3) 
  sb4_fold <- reIm2complex_sb(sb4) 
  
  fold_mat <-  sb1_fold
  fold_mat <- cbind(fold_mat, sb2_fold)
  fold_mat <- cbind(fold_mat, sb3_fold)
  fold_mat <- cbind(fold_mat, sb4_fold)
  

}



## Channel Reconstruction per 1 TTI

chan_reconstruction_row <- function(folded_mat, i, Nantennas=32){
  
  
  # Reconstruction of Channel H of all 1 TTI of MPCs
  
  # Inputs
  # folded_mat: matrix [MPCs x sb x TTI]
  
  # Outputs
  # sbfold: matrix [Nantennas  x Ncarriers * sb x TTI]
  
  # Intermediate Functions:
  # It calls 
  # Chan_reconstruction_row -> takes TTI by TTI and outputs the Channel
  
  
  ## SB 1 - Channel reconstruction
  SB1_chan = (matrix(G0 %*% folded_mat[,1,i], nrow = Nantennas))
  
  ## SB 2 - Channel reconstruction
  SB2_chan = (matrix(G0 %*% folded_mat[,2,i], nrow = Nantennas))
  
  ## SB 3 - Channel reconstruction
  SB3_chan = (matrix(G0 %*% folded_mat[,3,i], nrow = Nantennas))
  
  ## SB 4 - Channel reconstruction
  SB4_chan = (matrix(G0 %*% folded_mat[,4,i], nrow = Nantennas))
  
  ## Concat 4 SB channels
  
  fold_chan <-  SB1_chan
  fold_chan <- cbind(fold_chan, SB2_chan)
  fold_chan <- cbind(fold_chan, SB3_chan)
  fold_chan <- cbind(fold_chan, SB4_chan)
  
  fold_chan
  }
  
## Channel reconstruction of all the TTIs of Validation Set in 2D

channel_reconstruction <- function(folded_mat, Ncarriers=408, Nantennas=32){
  
  channel_all_TTI <- matrix(nrow = dim(folded_mat)[3] * Nantennas*4, ncol = Ncarriers) 
    
  for (i in 1:dim(folded_mat)[3]){
    
    channel <- chan_reconstruction_row(folded_mat,i)
    
    channel_all_TTI[(1+Nantennas*4*(i-1)):(Nantennas*4*i),] <- channel
  }
  
  channel_all_TTI
}

## Channel reconstruction of all the TTIs of Validation Set in 3D

channel_reconstruction_3dArray <- function(folded_mat, Ncarriers=408, Nantennas=32){
  
  # Reconstruction of Channel H of all TTIs in a folded Matrics of MPCsc
  
  # Inputs
  # folded_mat: matrix [MPCs x sb x TTI]
  
  # Outputs
  # sbfold: matrix [Nantennas x Ncarriers * sb x TTI]
  
  # Intermediate Functions:
  # It calls 
    # Chan_reconstruction_row -> takes TTI by TTI and outputs the Channel
  
  
  channel_all_TTI <- array(dim = c(Nantennas, Ncarriers*4, dim(folded_mat)[3]))
  
  for (i in 1:dim(folded_mat)[3]){
    
    channel <- chan_reconstruction_row(folded_mat,i)
    
    channel_all_TTI[,,i] <- channel
  }
  
  channel_all_TTI
}

## Divergence Metrics

mse_TTIs <- function(channel_pred, channel_truth){
  
  mse_vec <- rep(NA, dim(channel_pred)[3])
  
  for (i in 1:dim(channel_pred)[3]){
    
    ###   How do we take difference here ??
    
    mse_vec[i] <- mean((c(Mod(channel_truth[,,i])) - c(Mod(channel_pred[,,i])))^2)
    
  }
  
  mse_vec
}

#############################################################################
mse_TTIs_per_band <- function(channel_pred, channel_truth, Ncarriers=408){
  
  mse_vec <- matrix(nrow = 4 , ncol = dim(channel_pred)[3] )
  
  Sb_ind_1 = 1:Ncarriers
  Sb_ind_2 = (Ncarriers+1):(Ncarriers*2)
  Sb_ind_3 = (Ncarriers*2+1):(Ncarriers*3)
  Sb_ind_4 = (Ncarriers*3+1):(Ncarriers*4)
  
  for (i in 1:dim(channel_pred)[3]){
    
    ###   How do we take difference here ??
    
    mse_vec[ 1, i] <- norm(c(c(channel_truth[, Sb_ind_1, i]) - c(channel_pred[, Sb_ind_1, i])), type="2" )/norm(c(channel_truth[, Sb_ind_1,i]), type="2")
     
    mse_vec[ 2, i] <- norm(c(c(channel_truth[, Sb_ind_2, i]) - c(channel_pred[, Sb_ind_2, i])), type="2" )/norm(c(channel_truth[, Sb_ind_2,i]), type="2")
     
    mse_vec[ 3, i] <- norm(c(c(channel_truth[, Sb_ind_3, i]) - c(channel_pred[, Sb_ind_3, i])), type="2" )/norm(c(channel_truth[, Sb_ind_3,i]), type="2")
     
    mse_vec[ 4, i] <- norm(c(c(channel_truth[, Sb_ind_4, i]) - c(channel_pred[, Sb_ind_4, i])), type="2" )/norm(c(channel_truth[, Sb_ind_4,i]), type="2")
     
    
  }
  
  mse_vec
}


mse_TTIs_per_band_2 <- function(channel_pred, channel_truth, Ncarriers=408){
  
  mse_vec <- matrix(nrow = 4 , ncol = dim(channel_pred)[3] )
  
  Sb_ind_1 = 1:Ncarriers
  Sb_ind_2 = (Ncarriers+1):(Ncarriers*2)
  Sb_ind_3 = (Ncarriers*2+1):(Ncarriers*3)
  Sb_ind_4 = (Ncarriers*3+1):(Ncarriers*4)
  
  for (i in 1:dim(channel_pred)[3]){
    
    ###   How do we take difference here ??
    
    mse_vec[ 1, i] <- norm(((channel_truth[, Sb_ind_1, i]) - (channel_pred[, Sb_ind_1, i])), type="2" )/norm((channel_truth[, Sb_ind_1,i]), type="2")
    
    mse_vec[ 2, i] <- norm(((channel_truth[, Sb_ind_2, i]) - (channel_pred[, Sb_ind_2, i])), type="2" )/norm((channel_truth[, Sb_ind_2,i]), type="2")
    
    mse_vec[ 3, i] <- norm(((channel_truth[, Sb_ind_3, i]) - (channel_pred[, Sb_ind_3, i])), type="2" )/norm((channel_truth[, Sb_ind_3,i]), type="2")
    
    mse_vec[ 4, i] <- norm(((channel_truth[, Sb_ind_4, i]) - (channel_pred[, Sb_ind_4, i])), type="2" )/norm((channel_truth[, Sb_ind_4,i]), type="2")
    
    
  }
  
  mse_vec
}


#############################################################################
plot_MSE_by_band <- function(mse_mat, title, ylim=c(0,0.5)){
  
  sb <- dim(mse_mat)[1]
  
  plot(mse_mat[1,], main = title, ylim=ylim)
  
  for(i in 2:sb){
    
    points(mse_mat[i,], col = i)
    
  }
  
  legend(x = "bottom", legend=c("1st MPC SB1", "1st MPC SB2", "1st MPC SB3", "1st MPC SB4"),
         col=c( 1,2,3,4), lty = c(1,1), lwd = c(1,1))
  }

  


