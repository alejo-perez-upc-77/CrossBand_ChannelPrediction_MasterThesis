
library(R.matlab)
library(BiocManager)
library(pcaMethods)
library(MetabolAnalyze)
library("plot3D")
library(MASS)
library(abind)
library(plotly)
library(rgl)

source("../../Preprocessing_Utilities.R")
source("../../PPCA_Utilities.R")
source("../../FoldMatrix.R")


## Hardcoded Parameters

Ncarriers <- 408
Nantennas <- 32

#################################################
######## Channel Matrix Reconstruction   ########
#################################################

# We would Use Weight matrix here

data <- readMat("../../../Data/data_UE1_600/quadriga_ue1_snr10_avg10TTI_total600TTI_evalinterval1TTI_seed123.mat") 
  
## Weight Matrix Parameters of Transformation 

G0 <- data[1]$param[[1]]

mat <- data$param[[3]]


## Validation Prediction Set

pred_Pol_1 <- as.matrix(read.csv("./prediction.csv", row.names = 1))
pred_Pol_2 <- as.matrix(read.csv("./prediction_2ndpol.csv", row.names = 1))

folded_pred_pol1 <- fold_matrix(pred_Pol_1)

folded_pred_pol2 <- fold_matrix(pred_Pol_2)


## Folding Background Truth

valid_Set_TTI = dim(folded_pred_pol1)[3]

True_Pol_1 <- mat[1:(dim(mat)[1]/2), , (dim(mat)[3] - valid_Set_TTI + 1):(dim(mat)[3])]

True_Pol_2 <- mat[(dim(mat)[1]/2+1):dim(mat)[1], , (dim(mat)[3] - valid_Set_TTI + 1):(dim(mat)[3])]


################################################################################
#########  From here On a fix is needed. Generalize implementation 
##################           UE0 Doesn't Have       ############################
################################################################################



## Reconstruct the Predicted Matrix

Channel_total_POl_1 <- channel_reconstruction_3dArray(folded_pred_pol1)
Channel_total_POl_2 <- channel_reconstruction_3dArray(folded_pred_pol2)


Channel_True_Pol_1 <- channel_reconstruction_3dArray(True_Pol_1)
Channel_True_Pol_2 <- channel_reconstruction_3dArray(True_Pol_2)

## MSE Metric over TTIs


MSE_TTI <- mse_TTIs(Channel_total_POl_1, Channel_True_Pol_1)

mse_TTIs_per_band <- mse_TTIs_per_band(Channel_total_POl_1, Channel_True_Pol_1)


mse_TTIs_per_band[which(mse_TTIs_per_band<0.001)]<-NaN





plot_MSE_by_band(mse_TTIs_per_band)
  
#write.csv(folded_mat,'prediction_valid fold.csv')

plot_ly(z = mse_TTIs_per_band, type = "surface")


NTTI=dim(mse_TTIs_per_band)[2]

plot3d( 
  x=1:4, y=1:NTTI, z=mse_TTIs_per_band, 
  
  xlab="Subbands", ylab="TTIs", zlab="MSE")
