clear
close all

addpath('Toolbox/');

dataset = 'Ball';
do_write_obj = 0;	% Set to 0 to de-active OBJ writing, to 1 to activate it

%%% Set dataset
% File containing the initial lighting and intrinsics
calib_file = sprintf('Datasets/%s/calib.mat',dataset);
% Folder containing the photometric stereo images
images_folder = sprintf('Datasets/%s',dataset);

%%% Read calibration files
disp('Reading calibration files');
load(calib_file);
if(K(1,3)==0) K = K';end % Matlab uses a transposed definition of intrinsics matrix
% Store in a compact structure
calib.S = S; clear S; % Lighting (nimgs x 3 x nchannels)
calib.K = K; clear K; % Intrinsics (3 x 3)

%%% Read dataset
disp('Reading photometric stereo images');
% Get number of images from light calibration
nimgs = size(calib.S,1);
% Read ambient light image
Iamb = double(imread(sprintf('%s/ambient.png',images_folder)));
[nrows,ncols,nchannels] = size(Iamb);
% Read mask
mask = double(imread(sprintf('%s/mask.png',images_folder)));
mask = mask(:,:,1)>0;
% Read each image and substract ambient light
I = zeros(nrows,ncols,nchannels,nimgs);
for i = 1:nimgs
	Ii = double(imread(sprintf('%s/I%d.png',images_folder,i)));
	Ii = Ii-Iamb;
	I(:,:,:,i) = max(0,Ii);
end
clear i Ii Iamb light_file images_folder nrows ncols nchannels
% Store in a compact structure
data.I = I; clear I; % Images (nrows x ncols x nchannels x nimgs)
data.mask = mask; clear mask; % Images (nrows x ncols x nchannels x nimgs)

%%% Convert images to graylevel for this demo
data.I = squeeze(mean(data.I,3));
calib.S = mean(calib.S,3);

%%% Show the input images
disp('Displaying data');
figure(1001)
subplot(3,1,1)
imshow(uint8(255*data.I(:,:,1)./max(data.I(:))))
title('$$I^1$$','Interpreter','Latex','Fontsize',18);
subplot(3,1,2)
imshow(uint8(255*data.I(:,:,2)./max(data.I(:))))
title('$$I^2$$','Interpreter','Latex','Fontsize',18);
subplot(3,1,3)
imshow(uint8(255*data.mask))
title('$$\Omega$$','Interpreter','Latex','Fontsize',18);
drawnow

%%% Set parameters
disp('Setting parameters');
params.precond = 'cmg'; % Use multigrid preconditioner
params.z0 = 2000*double(data.mask); % Initial depth map: a plane at 2000mm from camera
params.estimator = 'Cauchy'; % Robust Cauchy M-estimator
params.ratio = 4; % Downsample by a factor of 4 for faster results 
params.self_shadows = 1; % Explicitly take into account self-shadows
params.display = 1; % Display result at each iteration
params.semi_calibrated = 0; % Do not refine intensities
params.uncalibrated = 1; % Refine whole lighting vectors
params.tol_pcg = 1e-9; % Stopping criterion for the inner CG iterations


%%% Solve photometric stereo
disp('Solving photometric stereo');
[XYZ,N,rho,Phi,S,mask,tab_nrj] = robust_ps(data,calib,params);
% Scale albedo for visualization
rho = rho./max(rho(:));
rho = uint8(255*rho);

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
