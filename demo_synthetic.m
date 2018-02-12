clear
close all

% Add the reconstruction functions
addpath(genpath('Toolbox/'));
% Set random seed to have comparable results in each testa
rng(31415926)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHETIC DATASET LOADING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load a synthetic dataset containing: 
%   I: nrows x ncols x nimgs (stack of graylevel images)
%   S: nimgs x 3 (light directions multiplied by intensities)
%   K: 3 x 3 (intrinsics matrix)
%   mask: nrows x ncols (binary mask of the object)
%   z: nrows x ncols (ground truth depth)
%   N: nrows x ncols x 3 (ground truth normals)
%   rho: nrows x ncols (ground truth albedo)
load('Datasets/synthetic_dataset.mat');

% Store a bunch of ground truth data for later
I_GT = I; max_I = max(I(:)); % Ground truth image and its maximum value
N_GT = N; % Ground truth normals
Phi_GT = sqrt(sum(S.^2,2)); % Ground truth light intensities
S_GT = bsxfun(@rdivide,S,Phi_GT); % Ground truth light directions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD VARIOUS PERTURBATIONS TO THE DATASET FOR ROBUSTNESS EVALUATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Noise added to the calibrated lighting intensities and directions (Gaussian)
std_intensity_bias = 25; % Standard deviation of the Gausian noise added to the light intensities (in percents of the max source intensity)
std_direction_bias = 30; % Standard deviation of the Gausian noise added to the light directions (in degrees)

% Noise added to the images (Gaussian + S&P outliers)
std_noise = 5; % Standard deviation of the Gausian noise (in percents of the max brightness)
prct_SP = 2; % Number of pixels corrupted by salt and pepper noise (in percents)

% Additive offset added to the image (constant)
% NOTE: such an ambient light map can be handled either:
%   1) by calibrating it and providing it to the algorithm in data.ambient (RECOMMENDED, only requires taking one more shot in real-world applications)
%   2) by letting the algorithm automatically estimate it, setting params.do_ambient = 1 (UNSTABLE)
offset = 30; % Additive offset map, aka ambient light (in percents of the max brightness) -- this 

% Call a script to add the perturbations. See add_perturbations.m 
add_perturbations; % Add calibration error and image noise 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY SOME DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SHOW THE FIRST TWO IMAGES
figure(101)
subplot(3,1,1)
imshow(uint8(255*I(:,:,1)./max(I(:))))
title('$$I^1$$','Interpreter','Latex','Fontsize',18);
subplot(3,1,2)
imshow(uint8(255*I(:,:,2)./max(I(:))))
title('$$I^2$$','Interpreter','Latex','Fontsize',18);
subplot(3,1,3)
imshow(uint8(255*mask))
title('Mask','Interpreter','Latex','Fontsize',18);
drawnow

%%% Display the perturbation statistics
disp(' ');
disp(' ');
disp('===================');
disp(' ');
disp('=== Statistics of the noise on the input images ===');
disp(sprintf('= Gaussian noise std: %.2f (in percents of the maximum uncorrupted images)',std_noise));
disp(sprintf('= Additive offset: %.2f (in percents of the maximum uncorrupted images)',offset));
disp(sprintf('= Number of outliers: %.2f (in percents)',prct_SP));
disp(' ');
disp('=== Standard deviation of the Gaussian noise on the calibrated lighting ===');
disp(sprintf('= Intensities: %.2f (in percents of the maxium uncorrupted intensity)',std_intensity_bias));
disp(sprintf('= Directions: %.2f (in degrees)',std_direction_bias));
disp(' ');
disp('===================');
disp(' ');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE DATASET FOR 3D-RECONSTRUCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Store dataset
data.I = I; % Images (nrows x ncols x nimgs, REQUIRED)
data.ambient = offset_map; % Ambient light  (nrows x nimgs - if not provided then it is assumed zero i.e. data.ambient = zeros(size(data.mask)))
data.mask = mask; % Binary mask (nrows x ncols - if not provided then the whole 2D domain is used i.e., data.mask = ones(size(I,1),size(I,2))))

%%% Store calibration data
calib.S = S; % Lighting directions (nimgs x 3, REQUIRED)
calib.Phi = Phi; % Lighting intensities (nimgs x 1, if not provided then calib.Phi = sqrt(sum((calib.S).^2,2)) is used)
calib.K = K; % Intrinsics (3 x 3, if not provided then orthographic camera is considered)

%%% Set parameters - See the script 'set_parameters.m' 
set_parameters; % Call a script to create a structure 'params' containing options for robust estimation, convergence criteria, initialization, etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALL MAIN ROUTINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');
disp(' ');
disp('===================');
disp(' ');
disp('=== Doing 3D-reconstruction ===');
[XYZ,N,rho,Phi,S,ambient,mask,tab_nrj] = robust_ps_V2(data,calib,params);
disp(' ');
disp('===================');
disp(' ');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Show the resulting shape
figure(1002);
surfl(-XYZ(:,:,1),-XYZ(:,:,2),-XYZ(:,:,3),[0 90]);
shading flat;
colormap gray
view(-220,30);
axis ij
axis image
title('Shape')

%%% Show the resulting albedo
figure(1003)
imagesc(rho,[0 max(I(:))])
colormap gray
colorbar
axis image
axis off
title('Albedo')

%%% Show the resulting normal map
figure(1004)
Ndisp = N;
Ndisp(:,:,3) = -Ndisp(:,:,3);
Ndisp = 0.5*(Ndisp+1);
imagesc(Ndisp)
axis image
axis off
title('Normals')
drawnow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVALUATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Compute angular error wrt ground truth normals (CAUTION: MEANINGFUL ONLY IF params.ratio = 1)
AE = rad2deg(real(acos(sum(N.*N_GT(1:params.ratio:end,1:params.ratio:end,:),3))));
mean_angular_error = mean(abs(AE(find(mask>0))));
median_angular_error = median(abs(AE(find(mask>0))));
std_angular_error = std(abs(AE(find(mask>0))));

%%% Show angular error map
figure(2000)
imagesc(AE,[0 45])
axis image
axis off
colorbar
title('Angular error map (in degrees');

disp(' ');
disp(' ');
disp('===================');
disp(' ');
disp('=== 3D-reconstruction errors ===');
disp(sprintf('= Mean absolute angular error on normals: %.2f (in degrees)',mean_angular_error));
disp(sprintf('= Median absolute angular error on normals: %.2f (in degrees)',median_angular_error));
disp(sprintf('= Std of absolute angular error on normals: %.2f (in degrees)',std_angular_error));
disp(' ');
disp('=== Mean absolute angular error on lighting directions ===');
disp(sprintf('= Initial: %.2f (in degrees)',mean(abs(rad2deg(acos(sum(calib.S.*S_GT,2)))))));
disp(sprintf('= Refined: %.2f (in degrees)',mean(abs(rad2deg(acos(sum(S.*S_GT,2)))))));
disp(' ');
disp('=== STD of the ratio of light intensities over ground truth ones  ===');
disp(sprintf('= Initial: %.2f (unit-less)',std(calib.Phi./Phi_GT)));
disp(sprintf('= Refined: %.2f (unit-less)',std(Phi./Phi_GT)));
disp(' ');
disp('===================');
disp(' ');

