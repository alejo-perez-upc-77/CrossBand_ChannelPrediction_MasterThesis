
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
####### Import Channel Background Truth   #######
#################################################

dataUE1 <- readMat("../../../Data/data_UE2_600/channel_ue2.mat") 

## UE1 Data  

True_chan_UE1_pol_1 <- dataUE1$channel.ue[1:Nantennas, , (10+1):dim(dataUE1$channel.ue)[3]]

True_chan_UE1_pol_2 <-  dataUE1$channel.ue[1:Nantennas, , (10+1):dim(dataUE1$channel.ue)[3]]


#################################################
######## Channel Matrix Reconstruction   ########
#################################################

# We would Use Weight matrix here

data <- readMat("../../../Data/data_UE2_600/quadriga_ue2_snr10_avg10TTI_total600TTI_evalinterval1TTI_seed123.mat") 
  
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

# Estimated_Pol_1 <- mat[1:(dim(mat)[1]/2), , (dim(mat)[3] - valid_Set_TTI + 1):(dim(mat)[3])]

# Estimated_Pol_2 <- mat[(dim(mat)[1]/2+1):dim(mat)[1], , (dim(mat)[3] - valid_Set_TTI + 1):(dim(mat)[3])]

Estimated_Pol_1 <- mat[1:(dim(mat)[1]/2), ,]

Estimated_Pol_2 <- mat[(dim(mat)[1]/2+1):dim(mat)[1], ,]

################################################################################
#########  From here On a fix is needed. Generalize implementation 
##################           UE0 Doesn't Have       ############################
################################################################################



## Reconstruct the Predicted Matrix

Channel_total_POl_1 <- channel_reconstruction_3dArray(folded_pred_pol1)
Channel_total_POl_2 <- channel_reconstruction_3dArray(folded_pred_pol2)


Channel_Estimated_Pol_1 <- channel_reconstruction_3dArray(Estimated_Pol_1)
Channel_Estimated_Pol_2 <- channel_reconstruction_3dArray(Estimated_Pol_2)

## MSE Metric over TTIs

################################################################################
#################### True Channel VS Estimated Channel #########################
################################################################################

MSE_TTI <- mse_TTIs(True_chan_UE1_pol_1, Channel_Estimated_Pol_1)

mse_TTIs_per_SB <- mse_TTIs_per_band(True_chan_UE1_pol_1,
                                       Channel_Estimated_Pol_1)

mse_TTIs_per_SB_2 <- mse_TTIs_per_band_2(True_chan_UE1_pol_1,
                                       Channel_Estimated_Pol_1)


mse_TTIs_per_band[which(mse_TTIs_per_band<0.0001)]<-NaN



######################################
######################################

plot_MSE_by_band(mse_TTIs_per_SB, title = 'NMSE of True Channel VS Estimated per band', ylim=c(0.065,0.115))
  
#write.csv(folded_mat,'prediction_valid fold.csv')

plot_ly(z = mse_TTIs_per_SB, type = "surface")


NTTI=dim(mse_TTIs_per_SB)[2]

plot3d( 
  x=1:4, y=1:NTTI, z=mse_TTIs_per_SB, 
  
  xlab="Subbands", ylab="TTIs", zlab="MSE")

######################################
######################################

plot_MSE_by_band(mse_TTIs_per_SB_2, title = 'NMSE of True Channel VS Estimated per band UE2', ylim=c(0.02,0.08))

#write.csv(folded_mat,'prediction_valid fold.csv')

plot_ly(z = mse_TTIs_per_SB_2, type = "surface")


NTTI=dim(mse_TTIs_per_SB_2)[2]

plot3d( 
  x=1:4, y=1:NTTI, z=mse_TTIs_per_SB_2, 
  
  xlab="Subbands", ylab="TTIs", zlab="MSE")

## 1 subpath





### Plot of one antenna


plot(Mod(True_chan_UE1_pol_1[1,1:100,1]), type='l', ylim=c(0,3))

points(Mod(Channel_Estimated_Pol_1[1,1:100,1]), type='l', col=2)







