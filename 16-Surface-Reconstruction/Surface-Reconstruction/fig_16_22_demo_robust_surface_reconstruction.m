%% Fig. 16.22 page 762
% Robust surface reconstruction
% Uwe Weidner's test image (ECCV, 1994)
%
% Wolfgang Förstner 07/14
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

addpath(genpath('../../General-Functions/'));
addpath('Functions')
plot_init;

close all

%% set parameters

init_rand = 6;             % init random generator - may be changed for other example

out_C = 0;    % 0 no covariance matrix as output
              % 1 covariance matrix as output

% output settings
print_type = 1;
plot_type  = 1;

%% plot settings
ss = plot_init;

%% initialize random number generation by fixed seed
init_rand_seed(init_rand);

%% ---------------------- LETS GO -----------------------------------------

disp('-----------------------------------------------------------')
disp('----- Fig. 16.22: Demo finding steps and break lines  -----')
disp('-----------------------------------------------------------')

%% generate dem point cloud
% [points,BB,delta_x,sigma_k,out_in,dem]=...
%     simulate_points_dem_6(70,10,10,1,0.00);

[points,BB,delta_x,sigma_k,out_in,dem]=...
    simulate_points_dem_6(70,10,10,1,0.00);

%% smooth dem non-robustly
type_robust = 0;
display('Nonrobust Smoothing ==============')
ds_nonrobust = smooth_dem_robust_bilinear...
    (points,BB,delta_x,sigma_k,out_C,type_robust,out_in,...
    print_type,plot_type);

%% smooth dem robustly
type_robust = 1;
display('Robust Smoothing =================')

[ds_robust,S,Sigma,Np,Nr,Mc,ver,A,w,w_f,W] = smooth_dem_robust_bilinear...
    (points,BB,delta_x,sigma_k,out_C,type_robust,out_in,...
    print_type,plot_type);
   

%% plot

% plot smoothed surface 
figure('name','Fig 16.22 Finding steps and break lines - smoothed dem','color','w',...
    'Position',[0.335*ss(1),0.52*ss(2),0.3*ss(1),0.4*ss(2)]);
plot_surface(ds_nonrobust,BB,delta_x,'plotfun',@mesh,'view',[-65,29]);
axis equal
title('smoothed dem','FontSize',16)

% plot residuals
figure('name','Fig 16.22 Finding steps and break lines - residuals of smoothed dem','color','w',...
    'Position',[0.335*ss(1),0.02*ss(2),0.3*ss(1),0.4*ss(2)]);
imagesc(ds_nonrobust-dem);
colormap(gray);
title('residuals of smoothed dem','FontSize',16)
axis equal; axis off;

% plot restored dem 
figure('name','Fig 16.22 Finding steps and break lines - restored dem','color','w',...
    'Position',[0.65*ss(1),0.52*ss(2),0.3*ss(1),0.4*ss(2)]);
plot_surface(ds_robust,BB,delta_x,'plotfun',@mesh,'view',[-65,29]);
axis equal
title('restored dem','FontSize',16)

% plot residuals
figure('name','Fig 16.22 Finding steps and break lines - residuals of restored dem','color','w',...
    'Position',[0.65*ss(1),0.02*ss(2),0.3*ss(1),0.4*ss(2)]);
imagesc(ds_robust-dem);
colormap(gray);
title('residuals of restored dem','FontSize',16)
axis equal; axis off;
