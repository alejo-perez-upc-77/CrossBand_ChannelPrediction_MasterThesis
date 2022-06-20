

unfold <- function(mat){
  
  # Unfolds Matrix "mat" with dimensions [p x sb x tti] into [tti x sb·p·2]. 
  # This 2 is stands for Re and Im Part Separated
  # p: Path
  # sb: Number of Freq Subband
  # tti: Time Transmitting Interval 
  
  # Inputs
  # mat: matrix [p x sb x tti]
  
  # Outputs
  # unfolded: matrix [tti x sb·p·2]
  
  p <- dim(mat)[1] # Paths per sb per Re or Im
  
  sb <- dim(mat)[2] 
  
  # empty matrix to be filled
  unfolded <- matrix(NA, dim(mat)[3], dim(mat)[1] * sb * 2)
  
  fold_iterate <- c(1:p, (2*p+1):(3*p), (4*p+1):(5*p), (6*p+1):(7*p))  
  # Iterates unfolded mat 
  
  fold_x <- 1
  
  for (i in 1:sb){ # Iterates sb (4)
    
    for  (j in 1:p){ # Iterates p
      
      subset <- mat[j,i,]
      
      unfolded[, fold_iterate[fold_x]]  <- Re(subset) 
      unfolded[, fold_iterate[fold_x] + p]  <- Im(subset) 
      
      fold_x <- fold_x + 1  
      
    }
  }
  return(unfolded)
}




mask_data <- function(mat,mask){
  
  # Masking of the data "mat" with a mask object given
  
  # Inputs
  # data: matrix [tti x sb·p·2]
  # mask: matrix [mask_length x sb·p·2]
  
  # Outputs
  # mat: matrix [tti x sb·p·2]
  
  bracket_mask <- nrow(mask)
  
  mat[(nrow(mat)-bracket_mask+1):nrow(mat), ] <-  mask * 
    mat[(nrow(mat)-bracket_mask+1):nrow(mat), ]
  mat
  
}

pilot_masking <- function(mat,perc){
  
  # It masks given a percentage of TTI beginning by last TTI. Periodically 
  # masking 75% of the paths except for pilot band
  # p_sb: Paths per sb per param x 2
  
  # Inputs
  # data: matrix [tti x sb·p·2]
  # perc: % of TTI masked 
  
  # Outputs
  # data_masked: matrix [tti x sb·p·2]
  
  p_sb <- ncol(mat)/4
  
  # find the subset of data rows that we have to mask
  
  masked_bracket_len <- floor(nrow(mat)*perc/4)*4
  
  # indexes for each masking pattern
  
  nomask_1 <- 1:p_sb
  nomask_2 <- (p_sb+1):(2*p_sb)
  nomask_3 <- (2*p_sb+1):(3*p_sb)
  nomask_4 <- (3*p_sb+1):ncol(mat)
  
  ones_row <- rep(1, ncol(mat))
  
  mask_1 <- ones_row
  mask_1[-nomask_1] <- NA
  
  mask_2 <- ones_row
  mask_2[-nomask_2] <- NA
  
  mask_3 <- ones_row
  mask_3[-nomask_3] <- NA
  
  mask_4 <- ones_row
  mask_4[-nomask_4] <- NA
  
  maskblock <- rbind(mask_1,  mask_3, mask_2, mask_4)
  
  multiplicator <- masked_bracket_len/4
  
  mask <- do.call(rbind, replicate(multiplicator, maskblock, simplify=FALSE)) 
  
  len_masked <<- nrow(mask)
  
  mask_data(mat, mask)
  
}


unfold_polar <- function(mat){
  
  # Unfolds Matrix "mat" with dimensions [p x sb x tti] into [tti x sb·p·2]. 
  # This 2 is stands for Re and Im Part Separated
  # p: Path
  # sb: Number of Freq Subband
  # tti: Time Transmitting Interval 
  
  # Inputs
  # mat: matrix [p x sb x tti]
  
  # Outputs
  # unfolded: matrix [tti x sb·p·2]
  
  p <- dim(mat)[1] # Paths per sb per Re or Im
  
  sb <- dim(mat)[2] 
  
  # empty matrix to be filled
  unfolded <- matrix(NA, dim(mat)[3], dim(mat)[1] * sb * 2)
  
  fold_iterate <- c(1:p, (2*p+1):(3*p), (4*p+1):(5*p), (6*p+1):(7*p))  
  # Iterates unfolded mat 
  
  fold_x <- 1
  
  for (i in 1:sb){ # Iterates sb (4)
    
    for  (j in 1:p){ # Iterates p
      
      subset <- mat[j,i,]
      
      unfolded[, fold_iterate[fold_x]]  <- Mod(subset) 
      unfolded[, fold_iterate[fold_x] + p]  <- Arg(subset) 
      
      fold_x <- fold_x + 1  
      
    }
  }
  return(unfolded)
}


