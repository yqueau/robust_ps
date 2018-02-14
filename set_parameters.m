%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS FOR AUTOMATIC REFINEMENT OF LIGHTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lighting intensities refinement ? 
params.semi_calibrated = 1; % Set to 1 to enable lighting intensity refinement (otherwise use the calibrated values provided in calib.Phi - default is 0)
% Lighting directions refinement ? 
params.uncalibrated = 1; % Set to 1 to enable lighting direction refinement (otherwise use the calibrated values provided in calib.S - default is 0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS FOR ROBUST ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Robust estimation ? 
params.estimator = 'Cauchy'; % can be 'LS', 'Cauchy', 'Lp', 'Welsh', 'Tukey' or 'GM' - default is least-squares ('LS')
% Robust estimator parameter
params.lambda = 0.15; % useless for 'LS', the 'p' value for Lp norm, otherwise see code for the advised value - e.g. 0.15 for Cauchy, 1 for Lp norm to use L1 norm optimization, etc.
% Self-shadows handling ? 
params.self_shadows = 1; % Set to 1 to explicitly take into account self-shadows (set to zero to use the linearized model without clamping shading to zero - default is 1)
% Explicit thresholding ? (may slightly improve results, if the threshold is set appropriately)
%~ params.used = (data.mask>0) & (data.I>0.01*max(data.I(:))) & (data.I<0.99*max(data.I(:))); % Mask of inliers (nrows x ncols x nimgs - if not provided then all the data are used i.e. params.used = data.mask


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET INITIALIZATION AND REGULARIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial depth map
params.z0 = 2000*ones(size((data.mask))); % Initial depth map: a plane at 2000mm from camera
% Depth prior ? 
params.prior = 1e-6; % Add a quadratic term (z-z0)^2, with weight params.prior. Large values are useful if z0 is already close to the solution, otherwise any small value is fine
% Non-convex stabilizer for the descent wrt the lighting
params.stabilizer = 0.1; % Add a penalty |s^(k)-s^(k-1)|^2, weighted by params.stabilizer, to prevent too large descent steps for the lighting
% Number of levels in the multi-scale pyramid
params.scales = 7; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STOPPING CRITERIA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stopping criterion 1
params.maxit = 100; % Stop if this number of iterations is reached
% Stopping criterion 2
params.tol = 1e-3; % Stop if the relative absolute variation of the photometric error is below this threshold
% Stopping criterion 3
params.tol_normals = 0.05; % Stop if the median angular difference between two consecutive estimates is below this threshold (in degrees)
% Preconditioner (CMG PRECONDITIONER CAN BE DOWNLOADED HERE: http://www.cs.cmu.edu/~jkoutis/cmg.html)
params.precond = 'jacobi'; % 'cmg' uses combinatorial multigrid (recommended), 'ichol' uses incomplete cholesky (not recommended)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS FOR THE INNER PCG ITERATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inner stopping criterion 1
params.maxit_pcg = 100; % Max number of inner PCG iterations
% Inner stopping criterion 2
params.tol_pcg = 1e-3; % Relative residual variation in inner PCG iterations  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data reduction
params.ratio = 1; % Downsample the data by a factor of this for faster results
% Display 
params.display = 1; % Display current result at each iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENTAL OPTIONS - USE WITH CAUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Albedo estimation ? WARNING: EXPERIMENTAL, USE WITH CAUTION (BETTER TO ALWAYS ESTIMATE IT)
params.do_depth = 1; % Set to 1 to enable depth estimation (otherwise assume depth fixed to params.z0 - default is 1)
% Albedo estimation ? WARNING: EXPERIMENTAL, USE WITH CAUTION (BETTER TO ALWAYS ESTIMATE IT)
params.do_albedo = 1; % Set to 1 to enable albedo estimation (otherwise assume constant albedo equal to 1 - default is 1)
% Ambient lighting estimation ? WARNING: EXPERIMENTAL, USE WITH CAUTION (BETTER TO CALIBRATE AMBIENT AND PROVIDE IT IN data.ambient)
params.do_ambient = 0; % Set to 1 to enable ambient lighting estimation (otherwise assume ambient light is given - default is 0)
