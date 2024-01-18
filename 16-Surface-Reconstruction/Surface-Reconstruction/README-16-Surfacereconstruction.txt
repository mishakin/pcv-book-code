README-file for 16-Surfacereconstruction
########################################

html-Documentation: call Doku-16-Surfacereconstruction/index.html

Main-files:
-----------

eq_15_120_theoretical_accuracy.m
fig_15_9b_demo_surface_reconstruction.m
fig_15_20_demo_surface_reconstruction.m
fig_15_20_quality_surface_reconstruction.m
fig_16_21_quality_surface_reconstruction.m
fig_15_22_Uwe_s_test.m
fig_15_23_aurelo_result.m



Procedures for reconstruction (in Folder Functions):
----------------------------------------------------
smooth_dem_robust_bilinear.m
smooth_dem_robust_bilinear_flat.m
estimate_dem_robust.m
estimate_dem_robust_bilinear.m
estimate_dem_robust_flat.m



Procedures for generating data (in Folder Functions):
-----------------------------------------------------
simulate_points_dem_0.m
simulate_points_dem_0_flat.m
simulate_points_dem_6.m
simulate_points_dem_10.m
simulate_points_dem_15_flat.m
simulate_points_dem_16_precision.m


Ancillary procedures (in Folder Functions):
-------------------------------------------
ply_read.m
ply_write.m
interpolate_bilinear.m


Data-file (in Folder Data):
---------------------------
fa2_aurelo_result_pyra0_ausgeschnitten-1.ply
fa2_aurelo_result_pyra0_ausgeschnitten-4.ply (larger set of points)


Sparse inverse (in Folder General-Functions):
---------------------------------------------
Sparse_inverse/sparseinv/sparseinv.m 

Routines for determining the sparse inverse of a matrix
(only those elements of the inverse which are non-zero before inversion)
Copyright (c) 2014, Tim Davis 

