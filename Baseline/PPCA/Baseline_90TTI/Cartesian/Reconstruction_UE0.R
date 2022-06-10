
library(R.matlab)
library(BiocManager)
library(pcaMethods)
library(MetabolAnalyze)
library("plot3D")
library(MASS)
library(abind)
library(plotly)
library(rgl)

source("../../../Preprocessing_Utilities.R")
source("../../../PPCA_Utilities.R")
source("../../../FoldMatrix.R")

#################################################
######## Channel Matrix Reconstruction   ########
#################################################

# We would Use Weight matrix here

#data <- readMat("../data_UE1_600/quadriga_ue1_snr10_avg10TTI_total600TTI_evalinterval1TTI_seed123.mat") 

## Weight Matrix Parameters of Transformation 

#G0 <- data[1]$param[[1]]



## Validation Prediction Set

pred_Pol_1 <- as.matrix(read.csv("./prediction.csv", row.names = 1))
pred_Pol_2 <- as.matrix(read.csv("./prediction_2ndpol.csv", row.names = 1))

folded_pred_pol1 <- fold_matrix(pred_Pol_1)

folded_pred_pol2 <- fold_matrix(pred_Pol_2)

## Folding Background Truth

data <- readMat("../../../../Data/data_90TTI/1ue_quandriga_90TTI_10dB_data.mat") 

mat <- data$amp

valid_Set_TTI = dim(folded_pred_pol1)[3]

True_Pol_1 <- mat[1:(dim(mat)[1]/2), , (dim(mat)[3] - valid_Set_TTI + 1):(dim(mat)[3])]

True_Pol_2 <- mat[(dim(mat)[1]/2+1):dim(mat)[1], , (dim(mat)[3] - valid_Set_TTI + 1):(dim(mat)[3])]


################################################################################
#########  From here On a fix is needed. Generalize implementation 
###################           UE0 Doesnt Have       ############################
################################################################################



## Reconstruct the Predicted Matrix

Channel_total_POl_1 <- channel_reconstruction_3dArray(folded_pred_pol1)
Channel_total_POl_2 <- channel_reconstruction_3dArray(folded_pred_pol2)


Channel_reconstruction_3dArray <- channel_reconstruction_3dArray(folded_mat_truth)

## MSE Metric over TTIs

mse_TTIs <- function(channel_pred, channel_truth){
  
  mse_vec <- rep(NA, dim(folded_mat)[3])
  
  for (i in 1:dim(folded_mat)[3]){
    
    mse_vec[i] <- mean((c(Mod(Channel_total[,,i])) - c(Mod(Channel_reconstruction_3dArray[,,i])))^2)
    
  }
  
  mse_vec
}

MSE_TTI <- mse_TTIs(Channel_total, Channel_reconstruction_3dArray)

#write.csv(folded_mat,'prediction_valid fold.csv')