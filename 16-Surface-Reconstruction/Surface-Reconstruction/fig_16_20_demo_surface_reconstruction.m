%% Fig. 16.20 page 761
% Surface reconstruction: test dem smoothing, plots of sparse matrices A, N etc.
%
% Wolfgang Förstner 07/14
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

addpath(genpath('../../General-Functions/'));
addpath('Functions')

close all

%% set parameters

% params for generating the data
init_rand = 6;            % seed for random generator - may be changed for other example

% smoothing
type_robust = 0;
%             0 not robust
%             1 only dem
%             2 only points
%             3 both
out_C = 1;    % 0 no covariance matrix as output
              % 1 covariance matrix as output

% output settings
print_type = 0;
plot_type  = 0;

%% initialize random number generation by fixed seed
init_rand_seed(init_rand);

%% ---------------------- LETS GO -----------------------------------------

disp('----- Fig. 16.20: Demo surface interpolation -----')

%% generate dem point cloud
[points,BB,delta_x,sigma_k,sigma_h,out_in,dem] = simulate_points_dem_0;

%% interpolate

[ds,S,Sigma,Np,Nr,Mc,ver,A,w,w_f,W] = smooth_dem_robust_bilinear...
    (points,BB,delta_x,sigma_k,out_C,type_robust,out_in,...
    print_type,plot_type);


%% plot

% plot settings
ss = plot_init;

% plot surface and points
figure('name','Fig 16.20 Demo surface interpolation','color','w',...
    'Position',[0.1*ss(1),0.2*ss(2),0.7*ss(1),0.6*ss(2)]);

subplot(1,2,1);hold on
plot_surface(ds,BB,delta_x,'alpha',0.3,'view',[-29,65],'colormap','cool');
for n=1:Np
    plot3(points(n,1),points(n,2),points(n,3),'.k','MarkerSize',15)
end
axis equal
xlabel('$x$');ylabel('$y$');zlabel('$z$');
title('Fitted surface','FontSize',16)

% plot standarddeviation
subplot(1,2,2);hold on
plot_surface(sqrt(S),BB,delta_x,'alpha',0.3,'view',[-29,65],'colormap','cool');
z=zeros(Np,1);
for n=1:Np
    z(n)=interpolate_bilinear(S,points(n,1),points(n,2),delta_x,BB,sigma_h);
    plot3(points(n,1),points(n,2),sqrt(z(n)),'.k','MarkerSize',15) 
end
plot3(0,0,0,'.k')
axis equal
xlabel('$x$');ylabel('$y$');zlabel('$\sigma_z$');
title('Standard deviations $\sigma_z$ of the estimated grid points','FontSize',16)



