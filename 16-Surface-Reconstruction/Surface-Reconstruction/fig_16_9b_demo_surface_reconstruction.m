%% Fig. 16.9 bottom page 740
% Reconstruct dem from 7 points with flattening and smoothing
%
% Wolfgang Förstner 01/15
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

addpath(genpath('../../General-Functions/'));
addpath('Functions')

close all

%% set parameters

% data generation
init_rand = 6;             % init random generator - may be changed for other example

% smoothing
type_robust = 0;
%             0 not robust
%             1 only dem
%             2 only points
%             3 both
out_C = 0;    % 0 no covariance matrix as output
              % 1 covariance matrix as output

% intermediate output
print_type = 1;
plot_type  = 0;


%% plot settings
ss = plot_init;

%% initialize random number generation by fixed seed
init_rand_seed(init_rand);


%% ---------------------- LETS GO -----------------------------------------

disp('----------------------------------------------------------------------------')
disp('----- demo reconstruct dem from 7 points with flattening and smoothing -----')
disp('----------------------------------------------------------------------------')

% grid size
dx=1/10;  % for final plot
%dx=1/4;
display(['Grid size = ',num2str(dx)])


disp(' ')
disp('... for final plot take grid size dx=1/10 (may take a minute)')
disp(' ')
 

%% generate dem point cloud
[points,BB,delta_x,sigma_k,sigma_s,dem,out_in] = simulate_points_dem_0_flat(dx);
   

%% interpolate flat
starttime_flattening = cputime;
ds_flat = smooth_dem_robust_bilinear_flat...
    (points,BB,delta_x,sigma_k,out_C,type_robust,...
    print_type,plot_type);

disp(['      time for solution flattening: ',num2str(cputime-starttime_flattening)])

%% interpolate smooth
starttime_smoothing = cputime;
ds_smooth = smooth_dem_robust_bilinear...
    (points,BB,delta_x,sigma_k,out_C,type_robust,...
    out_in,print_type,plot_type);    
    
disp(['      time for solution smoothing: ',num2str(cputime-starttime_smoothing)])


%% plot
figure('name','Fig 16.9 bottom: Best fitting 2D function','color','w',...
    'Position',[0.1*ss(1),0.2*ss(2),0.8*ss(1),0.5*ss(2)]);

subplot(1,2,1);hold on;
plot_surface(ds_flat,BB,delta_x,'EdgeColor','none','FaceLighting','gouraud',...
        'ColorFct','smoothtanh','shading','interp','alpha',0.3,'view',[-20,30]);
plot3(points(:,1),points(:,2),points(:,3),'.b','MarkerSize',15);
axis equal
title('flat dem with given points','FontSize',16)

subplot(1,2,2);hold on
plot_surface(ds_smooth,BB,delta_x,'EdgeColor','none','FaceLighting','gouraud',...
        'ColorFct','smoothtanh','shading','interp','alpha',0.3,'view',[-20,30]);   
plot3(points(:,1),points(:,2),points(:,3),'.b','MarkerSize',15)
axis equal
title('smooth dem with given points','FontSize',16)
  


