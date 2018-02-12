% 1) Add Gaussian noise (simulate sensor thermal noise) 
noise_map = std_noise*max_I*randn(size(I))/100;
I = I + noise_map; 

% 2) Add a systematic offset (simulate ambient light)
offset_map = offset*max_I*ones(size(mask))/100;
I = bsxfun(@plus,I,offset_map);

% 3) Salt and pepper noise (simulate cast shadows and specular spikes)
indices_px = reshape(1:length(I(:)),size(I,1)*size(I,2),size(I,3));
indices_px = indices_px(find(mask>0),:); % Indices of the masked pixels
nb_pixels_affected_by_sp = round(prct_SP*length(indices_px(:))/100); 
if(nb_pixels_affected_by_sp>0)
	pixels_affected_by_sp = indices_px(randperm(length(indices_px(:)),nb_pixels_affected_by_sp));
	I(pixels_affected_by_sp(1:round(0.5*nb_pixels_affected_by_sp))) = 0; % Cast shadows
	I(pixels_affected_by_sp(round(0.5*nb_pixels_affected_by_sp)+1:end)) = max_I+max(offset_map(:)); % Saturations
end

% 4) Bias the lighting intensities (simulate bad calibration)
intensity_bias = std_intensity_bias*randn(size(Phi_GT))*max(Phi_GT(:))/100; 
Phi = max(eps,Phi_GT+intensity_bias);

% 5) Bias the lighting directions (simulate bad calibration)
[th,ph,r] = cart2sph(S_GT(:,1),S_GT(:,2),S_GT(:,3)); 
th = th + deg2rad(std_direction_bias)*randn(size(th));
ph = min(0,max(-pi/2,ph + deg2rad(std_direction_bias)*randn(size(ph))));
[Sx,Sy,Sz] = sph2cart(th,ph,r); 
S = [Sx,Sy,Sz];
