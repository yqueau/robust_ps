# Codes_NonConvex_PS
Flexible Matlba codes for robust nonconvex photometric stereo, with various options including: 
- RGB or graylevel images
- orthographic or perspective camera
- Robustness to self-shadows, cast-shadows, highlights, etc. by: simple thresholding, inclusion of the self-shadows in the model, and robust M-estimation, etc.
- explicit self-shadows handling
- ability to refine lighting intensities (semi-calibrated PS), lighting directions, or both
- possibility to add a shape prior (e.g., rough Kinect or stereo estimate)
- ability to estimate an ambient light map

## Introduction

These Matlab codes implement the method for robust nonconvex photometric stereo described in [1]. Given a set of photometric stereo images of a still scene acquired from a fixed (pinhole or orthographic camera) camera, and (possibly inaccurate)  lighting parameters, this algorithm estimates depth, normals, albedo and, optionally, refined lighting intensities and directions. 

[1] "A Non-Convex Variational Approach to Photometric Stereo under Inaccurate Lighting", Yvain Quéau et al., Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (CVPR 2017).

Please cite the above work if using the provided codes for your own research. 

Author: Yvain Quéau, Technical University Munich, yvain.queau@tum.de


## Datasets and demo

A demo file on a synthetic dataset is provided
- `Datasets/synthetic_dataset.mat` contains a synthetic dataset 'synthetic_dataset.mat' containing a set of images of the 'Joyful Yell' 3D-model, along with ground truth shape, albedo, camera and lighting parameters. 
- `demo_synthetic.m` adds various types of corruptions to the images (noise, outliers, offset) and to the calibred lighting (intensities and directions), shows how to construct the arguments to the method, does the reconstruction and evaluates it. See the scripts `add_perturbations.m` and `set_parameters.m` for further details. Observe for instance what happens with this dataset if light directions are not refined, or if robust estimation is not carried out. 


## Usage

The main fuction is `Toolbox/robust_ps_V2.m` (see header for details).

Outputs:
- a gridded point cloud (the third dimension represents depth)
- normal map
- albedo
- lighting intensities
- lighting directions
- binary mask
- evolution of energy w.r.t. iterations 

Inputs: 

- a 'data' structure such that:
  * data.I contains the 3D or 4D image stack (**REQUIRED**)
  * data.mask contains a binary 2D mask
  * data.ambient contains an ambient light map

- a 'calib' structure such that:
  * calib.S contains the initial (possibly inaccurate) lighting vectors (**REQUIRED**)
  * calib.K contains the camera's intrinsics
  * calib.Phi contains the initial sources intensities (possibly inaccurate)

- a 'params' structure containing optional parameters. See below and the 'set_parameters.m' script for detailed advices. 

## Refining an inaccurate lighting calibration 

If the calibration of lighting directions and intensities is not perfect, it is possible to automatically refine them by turning on the following options:
  * params.semi_calibrated enables automatic intensity refinement
  * params.uncalibrated enables automatic lighting vectors refinement

## Enforcing robustness to shadows, highlights, etc.

Robustness can be improved by playing with the following parameters:
  * params.estimator sets the estimator. LS (least-squares) is a good initial choice, but robust M-estimators may be more accurate, though they require the parameter lambda to be set. The 'Cauchy' estimator is our favorite. 
  * params.self_shadows explicitly includes self-shadows in the model or not 
  * params.lambda sets the M-estimator parameter. For Cauchy estimator, lambda between 0.02 and 0.2 seems reasonable. Lp norm optimization is achieved by setting Lp as estimator, and p for lambda. See main function header for other possible choices
  * params.used might be used to provide the algorithm with a binary mask of the inliers (assuming you have a way to detect them)

## Providing a rough shape estimate 

If you have a good initial estimate of the shape (e.g., a Kinect prior), then you can provide it to the algorithm by using:
  * params.z0 is the initial depth. We advise to roughly (visually) estimate the distance from camera to object in mm, and set z0 to a constant matrix with this rough distance 
  * params.prior adds a quadratic regularization term enforcing the solution to stay close to params.z0, with weight params.smooth

## Controlling the non-convex descent

Various options can be set in order to control the behavior of the ARLS algorithm we put forward. Useful parameters are:
  * params.stabilizer adds a proximal gradient term in the descent wrt light (setting to high value is equivalent to not refining lighting, setting to low value may yield convergence issues due to the non-convex nature of the problem)
  * params.maxit sets the maximum number of outer iterations
  * params.tol sets the outer relative stopping criterion on the energy
  * params.tol_normals sets the stopping criterion on the median angular difference between two estimated normal maps, in degrees
  * params.maxit_pcg sets the max number of inner CG iterations within each global iteration
  * params.tol_pcg sets the inner CG iterations relative stopping criterion
  * params.precond sets the inner CG preconditioner. We strongly recommend to use `cmg` (see Dependencies), but if you want to stick to Matlab's builtin function, use `ichol`

## Fast debugging

For fast debugging or proof of concept, it may be useful to disable one or the other estimate, to reduce the size of the data or to display live surface, albedo, normals and energy:
  * params.ratio downsamples images by a factor of ratio
  * params.display displays the result at each iteration
  * params.do_depth indicates whether depth is refined or not (if not, depth is kept fixed to params.z0)
  * params.do_albedo indicates whether albedo is refined or not (if not, albedo is kept fixed to 1) 

## Ambient light estimation

In real-world scenarios, calibrating the ambient light is easy: just capture a series of images without any of the controlled lights on, average them and provide the result to data.ambient. If for some awkward reason you can't or don't want to do that, then ambient light can be estimated by setting the following option (WARNING: can be pretty unstable)
  * params.do_ambient indicates whether ambient light is refined or not (if not, ambient light is kept fixed to data.ambient)

## Dependencies

We strongly recommend to use the CMG preconditioner from Koutis et al., which can be downloaded here: 
http://www.cs.cmu.edu/~jkoutis/cmg.html

If CMG it is not installed, set the "precond" parameter to "ichol". Standard incomplete Cholesky will be used, which should prevent any error message in Matlab, but may also be super slow or even non-convergent. 





