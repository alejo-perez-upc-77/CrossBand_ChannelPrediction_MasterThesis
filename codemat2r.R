

W0 <- G0                                                                                                          
mySize <- size(channel[sub_band])                                                                                 
Hl0 <- reshape(channel[sub_band](1:ant_size/2,:),prod(mySize)/2,1)                                                
Hl1 <- reshape(channel[sub_band](ant_size/2+1:ant_size,:),prod(mySize)/2,1)                                       
pp_size <- length(pp)                                                                                             
pdp_power_path <- pp(1:pp_size,4)                                                                                 
diag_invSNR <- diag(1./pdp_power_path)                                                                            
ZF_W0 <- ((W0'*W0 + diag_invSNR) W0')                                                                             
a0 <- ZF_W0*Hl0                                                                                                   
a1 <- ZF_W0*Hl1                                                                                                  
UeResults.amp <- [a0;a1]                                                                                          
vTemp0 <- W0*a0                                                                                                   
vTemp1 <- W0*a1                                                                                                   
vTemp0 <- reshape(vTemp0,mySize(1)/2,mySize(2))                                                                   
vTemp1 <- reshape(vTemp1,mySize(1)/2,mySize(2))                                                                   
estimatedChJointUser <- [vTemp0;vTemp1]                                                                           
idealChannel <- squeeze(idealFreqChannel(:, (sub_band-1)*408+1:sub_band*408))                                 
mean_channelMSJointChEst <- mean(mean(abs(idealChannel-estimatedChJointUser).^2))/mean(mean(abs(idealChannel).^2))