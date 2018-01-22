% using simulated slice profiles, compare ghosting levels
% when we use the same RF pulse for each VFA excitation,
% versus our optimized pulses

N = 96; % image matrix size; 96 is divisible by 2,3,4

% load profiles from dzFleet
load mxyProfiles % contains Mxy and Mxy_sameRF
Nseg = size(Mxy,2); % number of EPI segments

% integrate slice profiles
Mxy_optRF_int = sum(Mxy);
Mxy_optRF_int = Mxy_optRF_int./mean(Mxy_optRF_int);
Mxy_sameRF_int = sum(Mxy_sameRF);
Mxy_sameRF_int = Mxy_sameRF_int./mean(Mxy_sameRF_int);

% get a phantom image
img = phantom(96);

% go to k-space and multiply each segment's signal by integrated profile;
% then go back to image space
img_sameRF = ft2(img);
img_optRF = ft2(img);
for ii = 1:Nseg
    img_sameRF(ii:Nseg:end,:) = img_sameRF(ii:Nseg:end,:)*Mxy_sameRF_int(ii);
    img_optRF(ii:Nseg:end,:) = img_optRF(ii:Nseg:end,:)*Mxy_optRF_int(ii);
end
img_sameRF = ift2(img_sameRF);
img_optRF = ift2(img_optRF);

% measure ghosting
ghostNRMSE_sameRF = norm(img_sameRF-img)...
    ./norm(img*mean(Mxy_sameRF_int));
ghostNRMSE_optRF = norm(img_optRF-img)...
    ./norm(img*mean(Mxy_optRF_int));

figure;
subplot(231)
imagesc(img);axis image;colormap gray,axis off;colorbar;title 'True Image'
subplot(232)
imagesc(abs(img_sameRF),[0 1]);axis image;colormap gray,axis off;colorbar
title 'Same RF'
subplot(233)
imagesc(abs(img_optRF),[0 1]);axis image;colormap gray,axis off;colorbar
title 'Optimized RF'
errMax = max([abs(img_sameRF(:)-img(:));...
    abs(img_optRF(:)-img(:))]);
subplot(235)
imagesc(abs(img_sameRF-img),[0 errMax]);axis image;
colormap gray,axis off;colorbar
title(sprintf('Same RF, Error. NRMSE = %0.1d%%',ghostNRMSE_sameRF*100));
subplot(236)
imagesc(abs(img_optRF-img),[0 errMax]);axis image;
colormap gray,axis off;colorbar
title(sprintf('Optimized RF, Error. NRMSE = %0.1d%%',ghostNRMSE_optRF*100));


