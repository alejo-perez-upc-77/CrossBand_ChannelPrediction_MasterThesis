library(R.matlab)
library(BiocManager)
library(pcaMethods)
library(MetabolAnalyze)
library("plot3D")
library(MASS)
library(abind)
library(plotly)
library(rgl)
source("./Baseline/Preprocessing_Utilities.R")

## Hardcoded Parameters

Ncarriers <- 408
Nantennas <- 32

## Import Channel data

dataUE1 <- readMat("./Data/data_UE1_600/channel_ue1.mat") 
dataUE2 <- readMat("./Data/data_UE2_600/channel_ue2.mat") 


## UE1 Data 

dataUE1_pol_1 <- dataUE1$channel.ue[1:Nantennas,, (10+1):dim(dataUE1$channel.ue)[3]]

dataUE1_pol_2 <-  dataUE2$channel.ue[1:Nantennas,, (10+1):dim(dataUE2$channel.ue)[3]]





