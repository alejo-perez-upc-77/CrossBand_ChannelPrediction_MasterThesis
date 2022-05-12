
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
  
  for (i in 1:p){
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

folded_mat <- fold_matrix(data_rec_postpro)

folded_mat_valid_unmask <- fold_matrix(valid_unmask)

folded_mat_train_reconstructed <- fold_matrix(train_reconstructed)


dim(folded_mat_train_reconstructed)

#write.csv(folded_mat,'prediction_valid fold.csv')