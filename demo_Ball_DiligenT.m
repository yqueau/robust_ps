clear
close all

addpath('Toolbox/');

dataset = 'Ball_DiliGenT';
do_write_obj = 0;	% Set to 0 to de-active OBJ writing, to 1 to activate it

%%% Set dataset
% File containing the initial lighting and intrinsics
calib_file = sprintf('Datasets/%s/calib.mat',dataset);
% Folder containing the photometric stereo images
images_folder = sprintf('Datasets/%s',dataset);

%%% Read calibration files
disp('Reading calibration files');
load(calib_file);
% Store in a compact structure
calib.S = S; clear S; % Lighting (nimgs x 3 x nchannels)
calib.K = K; clear K; % Intrinsics (3 x 3)

%%% Read dataset
disp('Reading photometric stereo images');
% Get number of images from light calibration
nimgs = size(calib.S,1);
% Read mask
mask = double(imread(sprintf('%s/mask.png',images_folder)));
mask = mask(:,:,1) >0;
[nrows,ncols] = size(mask);
nchannels = 3;
% Read each RGB image 
I = zeros(nrows,ncols,nchannels,nimgs);
for i = 1:nimgs
	Ii = double(imread(sprintf('%s/%03d.png',images_folder,i)));
	I(:,:,:,i) = max(0,Ii./(max(Ii(:))));
end
clear i Ii
% Store in a compact structure
data.I = I; clear I; % Images (nrows x ncols x nchannels x nimgs)
data.mask = mask; clear mask; % Images (nrows x ncols x nchannels x nimgs)

%%% Show the input images
disp('Displaying data');
figure(1001)
subplot(3,1,1)
imshow(uint8(255*data.I(:,:,:,1)./max(data.I(:))))
title('$$I^1$$','Interpreter','Latex','Fontsize',18);
subplot(3,1,2)
imshow(uint8(255*data.I(:,:,:,2)./max(data.I(:))))
title('$$I^2$$','Interpreter','Latex','Fontsize',18);
subplot(3,1,3)
imshow(uint8(255*data.mask))
title('$$\Omega$$','Interpreter','Latex','Fontsize',18);
drawnow


%%% Set parameters
disp('Setting parameters');
params.z0 = 700*double(data.mask); % Initial depth map
params.estimator = 'Lp'; % Robust Cauchy M-estimator
params.lambda = 0.2; % Parameter (here pseudo norm L^(0.2)
params.self_shadows = 1; % Explicitly take into account self-shadows
params.uncalibrated = 1; % Refine whole lighting vectors
params.display = 0; % Display result at each iteration
params.precond = 'ichol'; % Use multigrid preconditioner ('cmg' from Koutis is much better, but 'ichol' can be used with built-in Matlab functions)

%%% Convert images to graylevel
%data.I = squeeze(mean(data.I,3));
%calib.S = mean(calib.S,3);

%%% Solve photometric stereo
disp('Solving photometric stereo');
[XYZ,N,rho,Phi,S,mask,tab_nrj] = robust_ps(data,calib,params);
% Scale albedo for visualization
rho = rho./max(rho(:));

%%% Show the result
disp('Displaying results');
figure(1002);
surfl(-XYZ(:,:,1),-XYZ(:,:,2),-XYZ(:,:,3),[0 90]);
shading flat;
colormap gray
view(-220,30);
axis ij
axis image
title('Shape')

figure(1003)
imagesc(rho)
if(size(rho,3)==1)
	colormap gray
end
axis image
axis off
title('Albedo')

figure(1004)
Ndisp = N;
Ndisp(:,:,3) = -Ndisp(:,:,3);
Ndisp = 0.5*(Ndisp+1);
imagesc(Ndisp)
axis image
axis off
title('Normals')
drawnow

%%% Save to an obj file
if(do_write_obj)
	disp('Saving results');
	export_obj(XYZ,N,rho,mask,dataset);
end



% Compare with Ground truth
load(sprintf('%s/Ngt.mat',images_folder));
AE = rad2deg(real(acos(sum(N.*Ngt,3))));
figure(1005)
imagesc(AE,[0 30]);
colormap jet
colorbar
title('Agular error')

MAE = mean(AE(find(data.mask>0)));
disp(sprintf('MAE = %.3f deg',MAE));
