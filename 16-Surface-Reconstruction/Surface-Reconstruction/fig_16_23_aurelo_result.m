%% Fig. 16.23 page 764
% Reconstruction of the surface of a facade - dem from bundle adjustment result
% 
% Wolfgang Förstner 07/14
% last changes: Wolfgang Förstner 10/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de


close all
% clear all
clearvars

addpath(genpath('../../General-Functions/'));
addpath('Functions')
addpath('Data')

%% set parameters

init_rand = 6;             % init random generator - may be changed for other example

% resolution (number of surface grid elements at longer side)
resolution = 40;       % five times faster
resolution = 90;       % as in figure

% estimation
type_robust = 2;           % 0 not robust / 1 only dem / 2 only points / 3 both
out_C = 0;                 % 0 no covariance matrix as output / 1 covariance matrix as output

% output settings
print_type = 1;
plot_type  = 0;

%% plot settings
ss = plot_init;

%% initialize random number generation by fixed seed
init_rand_seed(init_rand);

%% ---------------------- LETS GO -----------------------------------------

disp('-----------------------------------------------------------------------------------------------------')
disp('----- Fig. 16.23: Reconstruction of the surface of a facade - dem from bundle adjustment result -----')
disp('-----------------------------------------------------------------------------------------------------')

display(['Number of grid points along longer side: ',num2str(resolution)])
display(' ')

%% generate dem point cloud
[points,BB,delta_x,sigma_k,tmp] = ...
    simulate_points_dem_10('fa2_aurelo_result_pyra0_ausgeschnitten-1.ply',resolution);

%% interpolate
starttime = cputime;
out_in = ones(size(points,1),1);
[ds,S,Sigma,Np,Nr,Mc,ver,A,w,w_f,W] = smooth_dem_robust_bilinear...
    (points,BB,delta_x,sigma_k,out_C,type_robust,out_in,print_type,plot_type);
disp(['      complete time for solution: ',num2str(cputime-starttime)])

%% plot

figure('name','Fig 16.23 dem from bundle adjustment result','color','w',...
    'Position',[0.5*ss(1),0.3*ss(2),0.4*ss(1),0.5*ss(2)]);
plot_surface(ds,BB,delta_x,'plotfun',@mesh,'view',[12,30]);
axis equal;axis off;
title('fitted dem - min curvature$^2$','FontSize',16)


