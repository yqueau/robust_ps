# Codes_NonConvex_PS
Code for robust nonconvex photometric stereo, under calibrated or inacurrate directional lighting.

## Introduction

These Matlab codes implement the method for robust nonconvex photometric stereo described in [1,2]. Given a set of photometric stereo images of a still scene acquired from a (pinhole or orthographic camera) camera, and (possibly inaccurate)  lightin parameters, this algorithm estimates depth, normals, albedo and, optionally, refined lightings. It can output a colored mesh in the .obj format. 

Features:
- Several datasets
- Graylevel or RGB-valued images
- Various robust estimation techniques (explicit self-shadows handling, least-squares or M-estimation, ...)
- Optional automatic refinement of lighting intensities (semi-calibrated setup), or of the whole lightings (inaccurate setup)
- Possibility to add control points
- Possibility to add box constraints on depth

[1] "A Non-Convex Variational Approach to Photometric Stereo under Inaccurate Lighting", Yvain Quéau et al., Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (CVPR 2017).

Please cite the above work if using the provided codes and/or datasets for your own research. 

Author: Yvain Quéau, Technical University Munich, yvain.queau@tum.de

Credits: Datasets provided by the Photosculptura company: https://www.photosculptura.nl/

## Datasets

- The `Datasets/` folder contains two datasets: a soccer ball painted with diffuse white painting, and a plastic doll. Each dataset contains 6 images, plus a calib.mat file containing roughly calibrated lightings and camera parameters. For the ball dataset, a mask is also provided. 

## Usage

The main fuction is `Toolbox/robust_ps.m` (see header for details).

Outputs:
- a gridded point cloud (the third dimension represents depth)
- normal map
- albedo
- lighting intensities
- lighting vectors
- binary mask
- evolution of energy w.r.t. iterations 

Inputs: 

- a data structure such that:
  * data.I contains the 3D or 4D image stack (**REQUIRED**)
  * data.mask contains a binary 2D mask
  * data.mask_control contains a binary 2D mask of control points

- a calib structure such that:
  * calib.S contains the initial lighting vectors (**REQUIRED**)
  * calib.K contains the camera's intrinsics
  * calib.Phi contains the initial sources intensities

- a param structure containing optional parameters. We strongly recommend to play a bit with the following parameters:
  * params.z0 is the initial depth. We advise to roughly (visually) estimate the distance from camera to object in mm, and set z0 to a constant matrix with this rough distance 
  * params.estimator sets the estimator. LS (least-squares) is a good initial choice, but robust M-estimators may be more accurate, though they require the parameter lambda to be set
  * params.lambda sets the M-estimator parameter. For Cauchy estimator, we use 0.1 in our datasets. L1 norm optimization is achieved by setting Lp as estimator, and 1 for lambda
  * params.self_shadows includes self-shadows in  the model or not
  * params.semi_calibrated enables automatic intensity refinement
  * params.uncalibrated enables automatic lighting vectors refinement
  * params.zmax sets max values for the depth, which may be useful to avoid artifacts
  * params.smooth adds a quadratic regularization term, with weight params.smooth

For fast debugging or proof of concept, it may be useful to reduce the size of the data, to limit the number of iterations, or to display live surface, albedo, normals and energy:
  * params.ratio downsamples images by a factor of ratio
  * params.maxit sets the maximum number of iterations
  * params.tol sets the relative stopping criterion on the energy
  * params.display displays the result at each iteration

The inner conjugate gradient iterations can be controlled by:
  * params.maxit_pcg sets the max number of CG iterations within each global iteration
  * params.tol_pcg sets its relative stopping criterion
  * params.precond sets the preconditioner. We strongly recommend to use `cmg` (see Dependencies), but if you want to stick to Matlab's builtin function, use `ichol`


## Demo

The following demo files are provided: 

- `demo_1_nonrobust_calibrated_ball.m` : demo of fast calibrated PS using graylevel-converted images, least-squares estimator, no self-shadows modeling, and image downsampling. 

- `demo_2_calibrated_color_robust_ps.m` : same, using robust estimation (Cauchy M-estimator + explicit self-shadows handling). Better results, but somewhat biased by the wrong calibration

- `demo_3_semicalibrated_color_robust_ps.m` : same, using automatic lighting refinement. Results are now quite good.  

- `demo_4_doll.m` : test on a more complicated dataset (no mask, specular, strong discontinuities...). Illustrates RGB images, box constraint, regularization, and control points (here, a boundary condition).
 

## Dependencies

We strongly recommend to use the CMG preconditioner from Koutis et al., which can be downloaded here: 
http://www.cs.cmu.edu/~jkoutis/cmg.html

If CMG it is not installed, set the "precond" parameter to "ichol". Standard incomplete Cholesky will be used, which should prevent any error message in Matlab, but may also be super slow or even non-convergent. 





