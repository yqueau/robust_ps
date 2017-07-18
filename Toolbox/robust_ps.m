function [XYZ,N,rho,Phi,S,mask,tab_nrj] = robust_ps(data,calib,params)
%NEAR_PS solves photometric stereo under nearby point light sources
%
%	=== USAGE ===	
%
%	XYZ = EXPORT_OBJ(DATA,CALIB) estimates an
%	NROWS x NCOLS x 3 gridded point cloud, using the images in DATA and
%	the calibration parameters in CALIB
%
%	XYZ = EXPORT_OBJ(DATA,CALIB,PARAMS) uses the algorithm parameters in
%	set in PARAMS
%
%	[XYZ,N] = EXPORT_OBJ(...) also provides an NROWS x NCOLS x 3 normal
%	map
%
%	[XYZ,N,rho] = EXPORT_OBJ(...) also provides an 
%	NROWS x NCOLS x NCHANNELS albedo map, where NCHANNELS = 1 if the
%	images are graylevel images, and NCHANNELS = 3 if they are RGB
%
%	[XYZ,N,rho,Phi] = EXPORT_OBJ(...) also provides an NIMGS x NCHANNELS
%	matrix containing lighting intensities
%
%	[XYZ,N,rho,Phi,S] = EXPORT_OBJ(...) also provides the lights
%
%	[XYZ,N,rho,Phi,S,mask] = EXPORT_OBJ(...) also provides the binary mask
%
%	[XYZ,N,rho,Phi,S,mask,tab_nrj] = EXPORT_OBJ(...) also provides a vector tab_nrj
%	containing the energy at each iteration
%
%	=== DATA ===
%	
%	- Required fields:
%		- DATA.I:	NROWS x NCOLS x NIMGS (graylevel images) or
%				NROWS x NCOLS x NCHANNELS x NIMGS (color images)
%	- Optional fields:
%		- DATA.mask:	M x N binary mask (default: ones(NROWS,NCOLS))
%		- DATA.mask_control:	M x N binary mask of the control points (default: ones(NROWS,NCOLS))
%
%	=== CALIB ===
%	- Required fields:
%		- CALIB.S:	NIMGS x 3 x NCHANNELS light sources
%
%	- Optional fields:
%		- CALIB.K:	3 x 3 intrinsics matrix
%				(default: [], i.e. orthographic)
%		- CALIB.Phi:	NIMGS x NCHANNELS intensities 
%				(default: ones(NIMGS,NCHANNELS)
%
%	=== PARAMETERS ===
%	- Optional fields:
%		- PARAMS.prior: prior parameter
%			(default: 0)
%		- PARAMS.smooth: smoothing parameter
%			(default: 0)
%		- PARAMS.semi_calibrated: binary variable indicating if
%		intensities are refined or not
%			(default: 0)
%		- PARAMS.uncalibrated: binary variable indicating if
%		lightings are refined or not
%			(default: 0)
%		- PARAMS.z0: NROWS x NCOLS initial depth map
%			(default: 700*ones(NROWS,NCOLS))
%		- PARAMS.zmax: NROWS x NCOLS max depth
%			(default: Inf*ones(NROWS,NCOLS))
%		- PARAMS.estimator: string indicating which estimator is used:
%		estimator = 'LS' uses least-squares
%		estimator = 'Cauchy' uses Cauchy's M-estimator
%		estimator = 'GM' uses Geman and McClure's M-estimator
%		estimator = 'Tukey' uses Tukey's biweight M-estimator
%		estimator = 'Welsh' uses Welsh's M-estimator
%			(default: 'LS')
%		- PARAMS.lambda: parameter of the robust estimator
%			(default: 0, since LS estimator requires no estimator)
%			(we recommend e.g. 0.1 for Cauchy)
%		- PARAMS.self_shadows: binary variable indicating if
%		self-shadows are included in the model or not	
%			(default: 1)
%		- PARAMS.maxit: max number of iterations
%			(default: 100)
%		- PARAMS.tol: relative stopping criterion on the energy
%			(default: 1e-3)
%		- PARAMS.ratio: integer downsampling factor 
%			(default: 1)
%		- PARAMS.display: binary variable indicating if current shape
%		and albedo are displayed at each iteration 
%			(default: 0)
%		- PARAMS.precond: preconditioner used
%		precond = 'cmg' uses Koutis conjugate multigrid preconditioner
%		precond = 'ichol' uses incomplete Cholesky 
%			(default: 'ichol', but 'cmg' is strongly recommended)
%		- PARAMS.tol_pcg: relative stopping criterion for inner CG iter
%			(default: 1e-4)
%		- PARAMS.maxit_pcg: max iterations for inner CG iter
%			(default: 50)
%
%	=== CREDITS ===
%
%	This is inspired by the robust PS method in :
%	[1]	"A Non-Convex Variational Approach to Photometric Stereo under Inaccurate Lighting"
%		Quéau et al.
%		Proc. CVPR 2017, Honolulu, USA 

	disp(' ');
	
	if(nargin<2)
		disp('ERROR: must provide data and calibration');
		return;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Check DATA
	% Check images
	if(~isfield(data,'I'))
		disp('ERROR: images not provided in DATA.I')
		return;
	end
	I = double(data.I); clear data.I;
	if(size(I,4)==1) % Graylevel images
		nchannels = 1;
		[nrows,ncols,nimgs] = size(I);
	else % Color images
		[nrows,ncols,nchannels,nimgs] = size(I);
	end
	% Check mask
	if(~isfield(data,'mask'))
		disp('WARNING: mask not provided in DATA.mask, using full 2D grid')
		data.mask = ones(nrows,ncols);
	end
	mask = double(data.mask)>0; clear data.mask;
	if(size(mask,1)~=nrows | size(mask,2)~=ncols | ndims(mask)>2)
		disp(sprintf('ERROR: mask should be %d x %d',nrows,ncols));
		return;
	end
	% Check mask
	if(~isfield(data,'mask_control'))
		disp('WARNING: mask_control not provided in DATA.mask_control, using full 2D grid')
		data.mask_control = zeros(nrows,ncols);
	end
	mask_control = double(data.mask_control)>0; clear data.mask_control;
	if(size(mask_control,1)~=nrows | size(mask_control,2)~=ncols | ndims(mask_control)>2)
		disp(sprintf('ERROR: mask_control should be %d x %d',nrows,ncols));
		return;
	end	
	clear data

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Check CALIB	
	% Check locations
	if(~isfield(calib,'S'))
		disp('ERROR: sources locations not provided in DATA.S')
		return;
	end
	S = double(calib.S); clear calib.S;
	if(size(S,1)~=nimgs | size(S,2)~= 3 | size(S,3)~= nchannels) 
		disp(sprintf('ERROR: S should be %d x 3 x %d',nimgs,nchannels));
		return;
	end
	% Check intrinsics
	if(~isfield(calib,'K')|size(calib.K,1)==0)
		disp('Warning: intrinsics not provided in CALIB.K - using orthographic projection')
		calib.K = [1,0,0.5*size(I,2);0,1,0.5*size(I,1);0,0,1];
		orthographic = 1;
	else
		orthographic = 0;
	end
	K = double(calib.K); clear calib.K;
	if(size(K,1)~=3 | size(K,2)~= 3 | ndims(K)>2) 
		disp('ERROR: K should be 3 x 3');
		return;
	end
	if(K(1,3)==0)
		K = K';
	end
	if(K(1,1)==0 | K(2,2)==0 | K(1,3)==0 | K(2,3)==0 | K(3,3)~=1 )
		disp('ERROR: intrinsics matrix not in a standard format');
		return;
	end
	% Check intensities
	if(~isfield(calib,'Phi'))
		disp('WARNING: intensities not provided in CALIB.Phi, using default values')
		calib.Phi = ones(nimgs,nchannels);
	end
	Phi = double(calib.Phi); clear calib.Phi;
	if(size(Phi,1)~=nimgs | size(Phi,2)~=nchannels | ndims(Phi)>2)
		disp(sprintf('ERROR: Phi should be %d x %d',nimgs,nchannels));
		return;
	end
	clear calib

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Check PARAMS
	if(~exist('params','var')|isempty(params)) params=[]; end;
	% Check prior
	if(~isfield(params,'prior'))
		disp('WARNING: prior not provided in PARAMS.prior, using default values')
		params.prior = 0e0;
	end
	prior = params.prior; clear params.prior;
	if(prior<0)
		disp(sprintf('ERROR: prior should be positive'));
	end
	% Check smooth
	if(~isfield(params,'smooth'))
		disp('WARNING: smooth not provided in PARAMS.smooth, using default values')
		params.smooth = 0e0;
	end
	smooth = params.smooth; clear params.smooth;
	if(smooth<0)
		disp(sprintf('ERROR: smooth should be positive'));
	end

	% Check initial depth
	if(~isfield(params,'z0'))
		disp('WARNING: initial depth not provided in PARAMS.z0, using default values')
		params.z0 = 700*ones(nrows,ncols);
	end
	z0 = params.z0; clear params.z0;
	if(size(z0,1)~=nrows | size(z0,2) ~= ncols | ndims(z0) > 2)
		disp(sprintf('ERROR: z0 should be %d x %d',nrows,ncols));
	end
	% Check min depth
	if(~isfield(params,'zmax'))
		disp('WARNING: max depth not provided in PARAMS.zmax, using default values')
		params.zmax = Inf*ones(size(mask));
	end
	zmax = params.zmax; clear params.zmax;
	if(size(zmax,1)~=nrows | size(zmax,2) ~= ncols | ndims(zmax) > 2)
		disp(sprintf('ERROR: zmax should be %d x %d',nrows,ncols));
		return
	end
	% Check semi_calibrated
	if(~isfield(params,'semi_calibrated'))
		disp('WARNING: semi_calibrated parameter not provided in PARAMS.semi_calibrated, using default values')
		params.semi_calibrated = 0;
	end
	semi_calibrated = params.semi_calibrated; clear params.semi_calibrated;
	if(semi_calibrated ~= 1 & semi_calibrated ~= 0)
		disp('ERROR: semi_calibrated should be binary');
	end
	% Check uncalibrated
	if(~isfield(params,'uncalibrated'))
		disp('WARNING: uncalibrated parameter not provided in PARAMS.uncalibrated, using default values')
		params.uncalibrated = 0;
	end
	uncalibrated = params.uncalibrated; clear params.uncalibrated;
	if(uncalibrated ~= 1 & uncalibrated ~= 0)
		disp('ERROR: uncalibrated should be binary');
	end	
	% Check estimator
	if(~isfield(params,'estimator'))
		disp('WARNING: estimator not provided in PARAMS.estimator, using default values')
		params.estimator = 'LS';
	end
	estimator = params.estimator; clear params.estimator;
	% Check estimator's parameter
	if(~isfield(params,'lambda'))
		disp('WARNING: lambda parameter not provided in PARAMS.lambda, using default values')
		if(strcmp(estimator,'LS'))
			params.lambda = 0;
		elseif(strcmp(estimator,'Cauchy'))
			params.lambda = 0.15;
		elseif(strcmp(estimator,'Lp'))
			params.lambda = 0.7;
		elseif(strcmp(estimator,'GM'))
			params.lambda = 0.4;
		elseif(strcmp(estimator,'Welsh'))
			params.lambda = 0.4;
		elseif(strcmp(estimator,'Tukey'))
			params.lambda = 0.9;
		end
	end
	lambda = params.lambda; clear params.lambda;
	if(~strcmp(estimator,'LS') & lambda == 0)
		disp('ERROR: a strictly positive estimator parameter must be set in PARAMS.lambda for robust M-estimation');
	end
	% Check self shadows
	if(~isfield(params,'self_shadows'))
		disp('WARNING: self_shadows parameter not provided in PARAMS.self_shadows, using default values')
		params.self_shadows = 1;
	end
	self_shadows = params.self_shadows; clear params.self_shadows;
	if(self_shadows ~= 1 & self_shadows ~= 0)
		disp('ERROR: self_shadows should be binary');
		return;
	end
	% Check max iterations
	if(~isfield(params,'maxit'))
		disp('WARNING: max number of iterations not set in PARAMS.maxit, using default values')
		params.maxit = 100;
	end
	maxit = params.maxit; clear params.maxit;
	if(maxit < 0)
		disp('ERROR: max number of iterations should be positive');
		return;
	end
	% Check tolerance
	if(~isfield(params,'tol'))
		disp('WARNING: tolerance not set in PARAMS.tol, using default values')
		params.tol = 1e-3;
	end
	tol = params.tol; clear params.tol;
	if(tol < 0)
		disp('ERROR: tolerance should be positive');
		return;
	end
	% Check ratio
	if(~isfield(params,'ratio'))
		disp('WARNING: downsampling factor not set in PARAMS.ratio, using default values')
		params.ratio = 1;
	end
	ratio = params.ratio; clear params.ratio;
	if(ratio < 0)
		disp('ERROR: ratio should be a positive integer');
		return;
	end
	if(ratio>1)
		K(1:2,:) = K(1:2,:)./ratio;
		mask = mask(1:ratio:end,1:ratio:end);
		mask_control = mask_control(1:ratio:end,1:ratio:end);
		z0 = z0(1:ratio:end,1:ratio:end);
		zmax = zmax(1:ratio:end,1:ratio:end);
		I = I(1:ratio:end,1:ratio:end,:,:);
		if(orthographic)
			K(1,1) = 1;
			K(2,2) = 1;
		end
		[nrows,ncols] = size(mask);
	end
	clear ratio
	% Check preconditioner
	if(~isfield(params,'precond'))
		disp('WARNING: preconditioner not provided in PARAMS.precond, using default values')
		params.precond = 'ichol';
	end
	precond = params.precond; clear params.precond;	
	if(~strcmp(precond,'cmg') & ~strcmp(precond,'ichol'))
		disp('ERROR: unknown preconditioner');
		return;
	end
	% Check display
	if(~isfield(params,'display'))
		disp('WARNING: display parameter not provided in PARAMS.display, using default values')
		params.display = 0;
	end
	display = params.display; clear params.display;
	if(display ~= 1 & display ~= 0)
		disp('ERROR: display should be binary');
		return;
	end
	if(display)
		hfig1 = figure();
		hfig2 = figure();
		hfig3 = figure();
		hfig4 = figure();
	end
	% Check PCG tolerance
	if(~isfield(params,'tol_pcg'))
		disp('WARNING: PCG tolerance not set in PARAMS.tol_pcg, using default values')
		params.tol_pcg = 1e-4;
	end
	tol_pcg = params.tol_pcg; clear params.tol_pcg;
	if(tol_pcg < 0)
		disp('ERROR: PCG tolerance should be positive');
		return;
	end
	% Check PCG max iterations
	if(~isfield(params,'maxit_pcg'))
		disp('WARNING: max number of iterations for PCG not set in PARAMS.maxit_pcg, using default values')
		params.maxit_pcg = 50;
	end
	maxit_pcg = params.maxit_pcg; clear params.maxit_pcg;
	if(maxit_pcg < 0)
		disp('ERROR: max number of iterations for PCG should be positive');
		return;
	end
	clear params
	disp(' ');

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Prepare data
	% Intrinsics
	fx = K(1,1);
	fy = K(2,2);
	x0 = K(1,3);
	y0 = K(2,3);
	clear K
	% Gradient operator on the mask
	G = make_gradient(mask);
	Dx =  G(1:2:end,:);
	Dy = G(2:2:end,:);
	Lap = Dx'*Dx+Dy'*Dy;
	clear G;
	% Estimator phi(x) and weight function w(x) = phi'(x)/x
	if(strcmp(estimator,'LS'))
		phi_fcn = @(x) x.^2;	
		w_fcn = @(x) 2*ones(size(x));	
	elseif(strcmp(estimator,'Cauchy'))
		phi_fcn = @(x) lambda^2*log(1+x.^2/lambda^2);	
	    w_fcn = @(x) 2*lambda^2./(x.^2+lambda^2);	
	elseif(strcmp(estimator,'Lp'))
		phi_fcn = @(x) (abs(x)).^lambda;
		thr_norm = 1e-2;
		w_fcn = @(x) lambda.*(max(thr_norm,abs(x))).^(lambda-2);
	elseif(strcmp(estimator,'GM'))
		phi_fcn = @(x) x.^2./(x.^2+lambda^2);
		w_fcn = @(x) 2*lambda^2./((x.^2+lambda^2).^2);	
	elseif(strcmp(estimator,'Welsh'))
		phi_fcn = @(x) lambda^2.*(1-exp(-x.^2./(lambda^2)));
		w_fcn = @(x) 2*exp(-x.^2./(lambda^2));
	elseif(strcmp(estimator,'Tukey'))
		phi_fcn = @(x) (abs(x)<=lambda).*(1-(1-x.^2./(lambda^2)).^3).*lambda^2+(abs(x)<lambda).*lambda^2;
		w_fcn = @(x) 6*(abs(x)<=lambda).*((1-x.^2./(lambda^2)).^2);
	end
	% Self shadow function psi(x) and derivative psi'(x)
	if(self_shadows)
		psi_fcn = @(x) max(x,0);     
		chi_fcn = @(x) double(x>=0); 
	else
		psi_fcn = @(x) x;   
		chi_fcn = @(x) 1*ones(size(x));	
	end
	% Divide images by lighting intensities and normalize to [0,1]
	if(nchannels == 1)
		for i = 1:nimgs
			I(:,:,i) = I(:,:,i)./Phi(i);
		end
	elseif(nchannels>1)
		if(size(Phi,2)==1)
			Phi = repmat(Phi,[1 nchannels]);
		end
		for i = 1:nimgs
			for ch = 1:nchannels
				I(:,:,ch,i) = I(:,:,ch,i)./Phi(i,ch);
			end
		end	
	end
	Phi(:) = 1;
	I = I./max(I(:));
	
	% Scaled pixel units
	[uu,vv] = meshgrid(1:ncols,1:nrows);
	u_tilde = (uu - x0); 
	v_tilde = (vv - y0);
	clear uu vv x0 y0
	% Control points
	imask = find(mask>0);
	imask_control = find(mask_control>0);
	ctrl = ismember(imask,imask_control);
	idx_control = find(ctrl>0);
	idx_ncontrol = find(ctrl == 0);
	npix = length(imask);
	npix_control = length(idx_control);
	clear imask_control mask_control ctrl
	% Some useful change
	if(nchannels>1)
		I = permute(I,[1 2 4 3]); 
	end
	% Some useful variables
	px_rep = repmat(u_tilde(imask),[1 nimgs]);
	py_rep = repmat(v_tilde(imask),[1 nimgs]);
	if(orthographic)
		px_rep(:) = 0;
		py_rep(:) = 0;
	end
	% Vectorize data
	I = reshape(I,nrows*ncols,nimgs,nchannels);
	I = I(imask,:,:);
	% Scale parameters
	if(~strcmp(estimator,'Lp'))
		lambda = lambda*median(abs(I(:)-median(I(:))));
	end
	smooth = nimgs*smooth;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Initialize variables
	disp(sprintf('Starting algorithm with %d pixels, %d images, %d channels, lambda = %.6f',npix,nimgs,nchannels,lambda));
	disp(' ')
	z = z0; 
	z(mask==0) = NaN;
	if(orthographic)
		z_tilde = z(imask);
		z0 = z0(imask);
		zmax = zmax(imask);
		XYZ = cat(3,u_tilde./fx,v_tilde./fy,z);
	else
		z0 = log(z0(imask));
		z_tilde = log(z(imask));
		zmax = log(zmax(imask));
		XYZ = cat(3,z.*u_tilde./fx,z.*v_tilde./fy,z);
	end
	rho = ones(nrows,ncols,nchannels);
	rho_tilde = ones(npix,nchannels); 
	zx = Dx*z_tilde;
	zy = Dy*z_tilde;
	Nx = zeros(nrows,ncols);
	Ny = zeros(nrows,ncols);
	Nz = zeros(nrows,ncols);
	Nx(imask) = fx*zx;
	Ny(imask) = fy*zy;
	if(orthographic)
		Nz(imask) = -1;
	else
		Nz(imask) = -u_tilde(imask).*zx-v_tilde(imask).*zy-1;
	end
	dz = sqrt(Nx.^2+Ny.^2+Nz.^2);
	N = cat(3,Nx./dz,Ny./dz,Nz./dz);
	tab_nrj = [];
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Initial energy
	shading_fcn = @(z,S) (spdiags(reshape(bsxfun(@minus,fx*transpose(S(:,1)),bsxfun(@times,px_rep,transpose(S(:,3)))),npix*nimgs,1),0,npix*nimgs,npix*nimgs)*repmat(Dx,[nimgs 1])+spdiags(reshape(bsxfun(@minus,fy*transpose(S(:,2)),bsxfun(@times,py_rep,transpose(S(:,3)))),npix*nimgs,1),0,npix*nimgs,npix*nimgs)*repmat(Dy,[nimgs 1]))*z-reshape(transpose(repmat(S(:,3),[1 npix])),npix*nimgs,1);
	r_fcn = @(rho,shadz,II) repmat(rho,[nimgs 1]).*psi_fcn(shadz)-II; % residual
	J_fcn = @(rho,shadz,II) sum(phi_fcn(r_fcn(rho,shadz,II))); % energy
	energy = 0;
	psi = zeros(npix,nimgs,nchannels);
	for ch = 1:nchannels
		psi(:,:,ch) = reshape(shading_fcn(z_tilde,S(:,:,ch)),npix,nimgs);
		psich = psi(:,:,ch);
		Ich = I(:,:,ch);
		energy = energy+J_fcn(rho_tilde(:,ch),psich(:),Ich(:));
	end
	clear Ich psich
	disp(sprintf('== it. 0 - energy : %.20f',energy));
	disp(' ');
	tab_nrj(1) = energy;
		

	% Display Initial result
	if(display)
		figure(hfig1);
		surfl(-XYZ(:,:,1),-XYZ(:,:,2),-XYZ(:,:,3),[0 90]);
		shading flat;
		colormap gray
		view(170,30);
		axis ij
		axis image
		title('Shape')

		figure(hfig2)
		semilogy(0:0,tab_nrj(1:1),'Linewidth',4)
		title('Energy')


		figure(hfig3)
		imagesc(rho)
		if(nchannels==1)
			colormap gray
			colorbar
		end
		axis image
		axis off
		title('Albedo')

		figure(hfig4)
		Ndisp = N;
		Ndisp(:,:,3) = -Ndisp(:,:,3);
		Ndisp = 0.5*(Ndisp+1);
		imagesc(Ndisp)
		axis image
		axis off
		title('Normals')
		drawnow
	end

	for it = 1:maxit
		w = zeros(npix,nimgs,nchannels);
		chi = chi_fcn(psi);

		% Intensities update
		if(semi_calibrated)
			% Phi update
			for ch = 1:nchannels
				Ich = I(:,:,ch);
				psich = psi(:,:,ch);
				w(:,:,ch) = w_fcn(reshape(r_fcn(rho_tilde(:,ch),psich(:),Ich(:)),npix,nimgs));
				rho_psi_chi = psi(:,:,ch).*chi(:,:,ch).*repmat(rho_tilde(:,ch),[1 nimgs]);
				Phi_rel(:,ch) = transpose(((sum(w(:,:,ch).*I(:,:,ch).*rho_psi_chi,1)))./(sum(w(:,:,ch).*(rho_psi_chi).^2,1)));
				I(:,:,ch) = bsxfun(@rdivide,I(:,:,ch),transpose(Phi_rel(:,ch)));
				Phi(:,ch) = Phi(:,ch).*Phi_rel(:,ch);
			end
			max_I = max(I(:)); 
			I = I./max_I;
			Phi = Phi./max_I;
		end
		clear Ich psich rho_psi_chi
		
		% Pseudo-albedo update
		for ch = 1:nchannels
			Ich = I(:,:,ch);
			psich = psi(:,:,ch);
			phi_chich = psi(:,:,ch).*chi(:,:,ch);
			w(:,:,ch) = w_fcn(reshape(r_fcn(rho_tilde(:,ch),psich(:),Ich(:)),npix,nimgs));
			denom = (sum(w(:,:,ch).*(phi_chich).^2,2));
			idx_ok = find(denom>0);
			rho_tilde(idx_ok,ch) = (sum(w(idx_ok,:,ch).*I(idx_ok,:,ch).*phi_chich(idx_ok,:),2))./denom(idx_ok);
		end
		clear psich phi_chich Ich

		% Log-depth update
		rho_rep = zeros(npix,nimgs*nchannels);
		for ch = 1:nchannels
			Ich = I(:,:,ch);
			psich = psi(:,:,ch);
			w(:,:,ch) = w_fcn(reshape(r_fcn(rho_tilde(:,ch),psich(:),Ich(:)),npix,nimgs));
			rho_rep(:,(ch-1)*nimgs+1:ch*nimgs) = repmat(rho_tilde(:,ch),[1 nimgs]);
		end
		
		D = reshape(chi,npix,nimgs*nchannels).*(rho_rep.^2).*reshape(w,npix,nimgs*nchannels);
		clear rho_rep

		A = [];
		rhs = [];
		for ch = 1:nchannels
			A =[A;(spdiags(reshape(bsxfun(@minus,fx*transpose(S(:,1,ch)),bsxfun(@times,px_rep,transpose(S(:,3,ch)))),npix*nimgs,1),0,npix*nimgs,npix*nimgs)*repmat(Dx,[nimgs 1])+spdiags(reshape(bsxfun(@minus,fy*transpose(S(:,2,ch)),bsxfun(@times,py_rep,transpose(S(:,3,ch)))),npix*nimgs,1),0,npix*nimgs,npix*nimgs)*repmat(Dy,[nimgs 1]))];
			rho_rep_ch = repmat(rho_tilde(:,ch),[1 nimgs]);
			%~ rhs_ch = chi(:,:,ch).*rho_rep_ch.*(rho_rep_ch.*(-psi(:,:,ch))+I(:,:,ch)).*w(:,:,ch);
			rhs_ch = chi(:,:,ch).*rho_rep_ch.*(bsxfun(@times,rho_rep_ch,transpose(S(:,3,ch)))+I(:,:,ch)).*w(:,:,ch);
			rhs = [rhs;rhs_ch(:)];
		end
		At = transpose(A);
		M = At*spdiags(D(:),0,nimgs*npix*nchannels,nimgs*npix*nchannels)*A;
		M = M+smooth*Lap+prior*speye(npix);		
		rhs = At*rhs(:)+prior*z0;
		clear A At rho_rep_ch rhs_ch
		
		% Set control points
		if(npix_control>0)
			%~ M = M(idx_ncontrol,idx_ncontrol);
			%~ rhs = rhs(idx_ncontrol);
			rhs = rhs-M(:,idx_control)*z_tilde(idx_control);
			rhs = rhs(idx_ncontrol);
			M = M(idx_ncontrol,idx_ncontrol);
		end	
		
		% Recompute the preconditioner every 5 iterations
		if(mod(it,5)==1)
			disp('Preconditioner recomputed')
			disp(' ');
			if(strcmp(precond,'cmg'))
				precond_L = cmg_sdd(M);
				precond_R = [];
			else
				precond_L = ichol(M);
				precond_R = precond_L';
			end
		end
		%~ z_tilde(idx_ncontrol) = z_tilde(idx_ncontrol)+ pcg(M,rhs,tol_pcg,maxit_pcg,precond_L,precond_R);
		z_tilde(idx_ncontrol) = pcg(M,rhs,tol_pcg,maxit_pcg,precond_L,precond_R,z_tilde(idx_ncontrol));

		% Enforce a box constraint
		z_tilde = min(z_tilde,zmax);

		% Auxiliary variables
		if(orthographic)
			z(imask) = z_tilde;
		else
			z(imask) = exp(z_tilde);
		end		
		zx = Dx*z_tilde;
		zy = Dy*z_tilde;
		Nx = zeros(nrows,ncols);
		Ny = zeros(nrows,ncols);
		Nz = zeros(nrows,ncols);
		Nx(imask) = fx*zx;
		Ny(imask) = fy*zy;
		if(orthographic)
			Nz(imask) = -1;
		else
			Nz(imask) = -u_tilde(imask).*zx-v_tilde(imask).*zy-1;
		end
		dz = sqrt(Nx.^2+Ny.^2+Nz.^2);
		N = cat(3,Nx./dz,Ny./dz,Nz./dz);
		for ch = 1:nchannels
			rhoch = rho(:,:,ch);
			rhoch(imask) = rho_tilde(:,ch).*dz(imask);
			max_val = median(rhoch(imask))+8*std(rhoch(imask));
			rhoch(rhoch>max_val) = max_val;
			rhoch(rhoch<0) = 0;
			rho(:,:,ch) = rhoch;
		end
		if(orthographic)
			XYZ = cat(3,u_tilde,v_tilde,z);
		else
			XYZ = cat(3,z.*u_tilde./fx,z.*v_tilde./fy,z);
		end
		
		%%% Update S
		if(uncalibrated & it>3)
			NNt_11 = Nx(imask).^2;
			NNt_12 = Nx(imask).*Ny(imask);
			NNt_13 = Nx(imask).*Nz(imask);
			NNt_22 = Ny(imask).^2;
			NNt_23 = Ny(imask).*Nz(imask);
			NNt_33 = Nz(imask).^2;
			
			for ch = 1:nchannels
				Ich = I(:,:,ch);
				psich = psi(:,:,ch);
				rho_rep = repmat(rho_tilde(:,ch),[1 nimgs]);
				w(:,:,ch) = w_fcn(reshape(r_fcn(rho_tilde(:,ch),psich(:),Ich(:)),npix,nimgs));
				clear psich
				% Aliases
				w_rho_chi = w(:,:,ch).*rho_rep.*chi(:,:,ch);
				w_rho2_chi2 = w_rho_chi.*rho_rep.*chi(:,:,ch);
				clear rho_rep
				
				% Construct second members
				w_rho_chi_I = w_rho_chi.*I(:,:,ch);
				clear w_rho_chi
				w_rho_chi_I_N = [sum(bsxfun(@times,w_rho_chi_I,Nx(imask)),1);
								sum(bsxfun(@times,w_rho_chi_I,Ny(imask)),1);...
								sum(bsxfun(@times,w_rho_chi_I,Nz(imask)),1)];
				sec_mb =  reshape(w_rho_chi_I_N,[3*nimgs,1]);
				clear w_rho_chi_I w_rho_chi_I_N
				
				% Construct matrices 
				mat_i_11 = sum(bsxfun(@times,w_rho2_chi2,NNt_11),1);
				mat_i_12 = sum(bsxfun(@times,w_rho2_chi2,NNt_12),1);
				mat_i_13 = sum(bsxfun(@times,w_rho2_chi2,NNt_13),1);
				mat_i_22 = sum(bsxfun(@times,w_rho2_chi2,NNt_22),1);
				mat_i_23 = sum(bsxfun(@times,w_rho2_chi2,NNt_23),1);
				mat_i_33 = sum(bsxfun(@times,w_rho2_chi2,NNt_33),1);
				clear w_rho2_chi2
				rows_mat_i = [1:3:3*nimgs-2 1:3:3*nimgs-2 1:3:3*nimgs-2 2:3:3*nimgs-1 2:3:3*nimgs-1 2:3:3*nimgs-1 3:3:3*nimgs 3:3:3*nimgs 3:3:3*nimgs];
				cols_mat_i = [1:3:3*nimgs-2 2:3:3*nimgs-1 3:3:3*nimgs 1:3:3*nimgs-2 2:3:3*nimgs-1 3:3:3*nimgs 1:3:3*nimgs-2 2:3:3*nimgs-1 3:3:3*nimgs];
				mat_i = sparse(rows_mat_i,cols_mat_i,[mat_i_11 mat_i_12 mat_i_13 mat_i_12 mat_i_22 mat_i_23 mat_i_13 mat_i_23 mat_i_33]); 
				
				S(:,:,ch) = transpose(reshape(mat_i\sec_mb,[3 nimgs]));
				psi(:,:,ch) = reshape(shading_fcn(z_tilde,S(:,:,ch)),npix,nimgs);
			end
			% update chi
			chi = chi_fcn(psi);
		end
		

		% Convergence test
		energy_new = 0;
		for ch = 1:nchannels
			psi(:,:,ch) = reshape(shading_fcn(z_tilde,S(:,:,ch)),npix,nimgs);
			Ich = I(:,:,ch);
			psich = psi(:,:,ch);
			energy_new = energy_new+J_fcn(rho_tilde(:,ch),psich(:),Ich(:));
		end
		relative_diff = abs(energy_new-energy)./energy_new;
		disp(sprintf('== it. %d - energy : %.6f - relative diff : %.6f',it,energy_new,relative_diff));
		energy = energy_new;
		tab_nrj(it+1) = energy;

		% Display current result
		if(display)
			figure(hfig1);
			surfl(-XYZ(:,:,1),-XYZ(:,:,2),-XYZ(:,:,3),[0 90]);
			shading flat;
			colormap gray
			view(170,30);
			axis ij
			axis image
			title('Shape')

			figure(hfig2)
			semilogy(0:it,tab_nrj(1:it+1),'Linewidth',4)
			title('Energy')

			figure(hfig3)
			rho_disp = rho./max(rho(:));
			imagesc(uint8(255*rho_disp))
			if(nchannels==1)
				colormap gray
				colorbar
			end
			axis image
			axis off
			title('Albedo')

			figure(hfig4)
			Ndisp = N;
			Ndisp(:,:,3) = -Ndisp(:,:,3);
			Ndisp = 0.5*(Ndisp+1);
			imagesc(Ndisp)
			axis image
			axis off
			title('Normals')
						
			drawnow	
		end
		disp(' ');
		if(relative_diff<tol & it > 5)
			break;
		end
	end
	if(display)
		close(hfig1);
		close(hfig2);
		close(hfig3);
		close(hfig4);
	end
	disp(' ');
end	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	

% Functions for computing the gradient operator on non-rectangular domains
function [M,Dyp,Dym,Dxp,Dxm] = make_gradient(mask)

	% Compute forward (Dxp and Dyp) and backward (Dxm and Dym) operators
	[Dyp,Dym,Dxp,Dxm,Sup,Sum,Svp,Svm,Omega,index_matrix,imask] = gradient_operators(mask);
	[nrows,ncols] = size(mask);

	% When there is no bottom neighbor, replace by backward (or by 0 if no top)
	Dy = Dyp;
	no_bottom = find(~Omega(:,:,1));
	no_bottom = nonzeros(index_matrix(no_bottom));
	Dy(no_bottom,:) = Dym(no_bottom,:);

	% Same for the x direction (right / left)
	Dx = Dxp;
	no_right = find(~Omega(:,:,3));
	no_right = nonzeros(index_matrix(no_right));
	Dx(no_right,:) = Dxm(no_right,:);

	M = sparse([],[],[],2*size(Dx,1),size(Dx,2),2*length(imask));
	M(1:2:end-1,:) = Dx;
	M(2:2:end,:) = Dy;
end

function [Dup,Dum,Dvp,Dvm,Sup,Sum,Svp,Svm,Omega,index_matrix,imask] = gradient_operators(mask)

	[nrows,ncols] = size(mask);
	Omega_padded = padarray(mask,[1 1],0);

	% Pixels who have bottom neighbor in mask
	Omega(:,:,1) = mask.*Omega_padded(3:end,2:end-1);
	% Pixels who have top neighbor in mask
	Omega(:,:,2) = mask.*Omega_padded(1:end-2,2:end-1);
	% Pixels who have right neighbor in mask
	Omega(:,:,3) = mask.*Omega_padded(2:end-1,3:end);
	% Pixels who have left neighbor in mask
	Omega(:,:,4) = mask.*Omega_padded(2:end-1,1:end-2);
	

	imask = find(mask>0);
	index_matrix = zeros(nrows,ncols);
	index_matrix(imask) = 1:length(imask);

	% Dv matrix
	% When there is a neighbor on the right : forward differences
	idx_c = find(Omega(:,:,3)>0);
	[xc,yc] = ind2sub(size(mask),idx_c);
	indices_centre = index_matrix(idx_c);	
	indices_right = index_matrix(sub2ind(size(mask),xc,yc+1));
	indices_right = indices_right(:);
	II = indices_centre;
	JJ = indices_right;
	KK = ones(length(indices_centre),1);
	II = [II;indices_centre];
	JJ = [JJ;indices_centre];
	KK = [KK;-ones(length(indices_centre),1)];
	
	Dvp = sparse(II,JJ,KK,length(imask),length(imask));
	Svp = speye(length(imask));
	Svp = Svp(index_matrix(idx_c),:);

	% When there is a neighbor on the left : backward differences
	idx_c = find(Omega(:,:,4)>0);
	[xc,yc] = ind2sub(size(mask),idx_c);
	indices_centre = index_matrix(idx_c);	
	indices_right = index_matrix(sub2ind(size(mask),xc,yc-1));
	indices_right = indices_right(:);
	II = [indices_centre];
	JJ = [indices_right];
	KK = [-ones(length(indices_centre),1)];
	II = [II;indices_centre];
	JJ = [JJ;indices_centre];
	KK = [KK;ones(length(indices_centre),1)];
	
	Dvm = sparse(II,JJ,KK,length(imask),length(imask));
	Svm = speye(length(imask));
	Svm = Svm(index_matrix(idx_c),:);



	% Du matrix
	% When there is a neighbor on the bottom : forward differences
	idx_c = find(Omega(:,:,1)>0);
	[xc,yc] = ind2sub(size(mask),idx_c);
	indices_centre = index_matrix(idx_c);	
	indices_right = index_matrix(sub2ind(size(mask),xc+1,yc));
	indices_right = indices_right(:);
	II = indices_centre;
	JJ = indices_right;
	KK = ones(length(indices_centre),1);
	II = [II;indices_centre];
	JJ = [JJ;indices_centre];
	KK = [KK;-ones(length(indices_centre),1)];
	
	Dup = sparse(II,JJ,KK,length(imask),length(imask));
	Sup = speye(length(imask));
	Sup = Sup(index_matrix(idx_c),:);

	% When there is a neighbor on the top : backward differences
	idx_c = find(Omega(:,:,2)>0);
	[xc,yc] = ind2sub(size(mask),idx_c);
	indices_centre = index_matrix(idx_c);	
	indices_right = index_matrix(sub2ind(size(mask),xc-1,yc));
	indices_right = indices_right(:);
	II = [indices_centre];
	JJ = [indices_right];
	KK = [-ones(length(indices_centre),1)];
	II = [II;indices_centre];
	JJ = [JJ;indices_centre];
	KK = [KK;ones(length(indices_centre),1)];
	
	Dum = sparse(II,JJ,KK,length(imask),length(imask));
	Sum = speye(length(imask));
	Sum = Sum(index_matrix(idx_c),:);	

end
