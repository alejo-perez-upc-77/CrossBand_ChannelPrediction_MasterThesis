W0 = G0; #param.G0
mySize = size(channel{sub_band}); # channel is the noise contaminated channel
Hl0 = reshape(channel{sub_band}(1:ant_size/2,:),prod(mySize)/2,1); # ant_size = 64, the 1st polarization
Hl1 = reshape(channel{sub_band}(ant_size/2+1:ant_size,:),prod(mySize)/2,1); # the 2nd polarization

pp_size = length(pp);
#ZF_W0 = ((W0'*W0)\W0');
#if ~SimPara.SRSChEst_clustering
pdp_power_path = pp(1:pp_size,4); # power of each peak of pdp (after denoise)
diag_invSNR = diag(1./pdp_power_path); # big path - high SNR, so "less noise" -- "noise level"
ZF_W0 = ((W0'*W0 + diag_invSNR)\W0'); # (zero forcing if without diag_invSNR) model the noise (sigma in gaussian)
a0 = ZF_W0*Hl0; # amplitudes of peaks in polarization 1
a1 = ZF_W0*Hl1; # amplitudes of peaks in polarization 2
UeResults.amp = [a0;a1]; ### stack the amplitudes of both polarizations

vTemp0 = W0*a0; #reconstruction from amplitude to channel 
vTemp1 = W0*a1;
vTemp0 = reshape(vTemp0,mySize(1)/2,mySize(2));
vTemp1 = reshape(vTemp1,mySize(1)/2,mySize(2));
estimatedChJointUser = [vTemp0;vTemp1];
idealChannel = squeeze(idealFreqChannel(:, (sub_band-1)*408+1:sub_band*408));
mean_channelMSJointChEst = mean(mean(abs(idealChannel-estimatedChJointUser).^2))/mean(mean(abs(idealChannel).^2));


