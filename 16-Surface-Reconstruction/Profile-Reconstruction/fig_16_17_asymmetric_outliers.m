%% Fig. 16.17 page 756
% Reconstruction of a profile with one-sided outliers
%
% Wolfgang Förstner 2014-08-06
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

close all
% clear all
clearvars

addpath(genpath('../../General-Functions/'));
addpath('Functions')


%% plot settings
ss = plot_init;

disp('----- eliminate asymmetric outliers -----')

%% set parameters

% params for generating the data
init_rand = 8;     % may be changed for other example
N = 200;           % number of points
dens = 0.8;        % percentage of observed points in (0,1]
sigma_e = 0.5;     % process noise
sigma_n = 0.5;     % observation noise
Poutlier = 0.4;    % fraction of outliers
Max_outlier = 25;  % max(outlier)
type_outlier = 1;  % asymmetric %type_outlier = 0; % symmetric

% params for estimation
type_robust  = [1,2,2,2]; %  
    % 1. [0 = L1, 1 = Kraus,
    % 2. g_factor,
    % 3. w_factor
    % 4. 0,1,2,3 = 00,01,10,11 robust for points and dem

print_type = 0;
plot_type  = 0;

Niter = 6;
if type_outlier ==0
    factor_sigma = 1;
else
    factor_sigma = 8;
end


%% initialize random number generation by fixed seed
init_rand_seed(init_rand);

%% generate profile
[x,y,out_in,select,xs,ys] = ...
    generate_observed_AR2(N,sigma_e,sigma_n,Poutlier,Max_outlier,type_outlier,dens);


%% reconstruct profile L1
type_robust(1) = 0;


xest = estimate_profile_robust...
    (N,select,ys,sigma_e/factor_sigma,sigma_n,Niter,type_outlier,...
    type_robust,print_type,plot_type);


%% reconstruct profile Kraus
type_robust(1) = 1;
factor_sigma= 1;

xest_kraus = estimate_profile_robust...
    (N,select,ys,sigma_e/factor_sigma,sigma_n,Niter,type_outlier,...
    type_robust,print_type,plot_type);


%% plot
figure('name','Fig. 16.17: Reconstruction of a profile with one-sided outliers',...
    'color','w', 'Position',[0.2*ss(1),0.1*ss(2),0.4*ss(1),0.8*ss(2)]);

subplot(2,1,1);hold on
plot(1:N,x,'--r','LineWidth',3)
plot(1:N,xest,'-k','LineWidth',2)
plot(select,ys,'.b','MarkerSize',20)
xlabel('$x$');ylabel('$z$');
title('Fig. 16.17a: Reconstruction with one-sided $L_{12}$-weight function');

subplot(2,1,2);hold on
plot(select,ys,'.b','MarkerSize',20)
plot(1:N,x,'--r','LineWidth',3)
plot(1:N,xest_kraus,'-k','LineWidth',2)
xlabel('$x$');ylabel('$z$');
title(['Fig. 16.17b: Reconstruction with the weight function of Kraus and Pfeifer,' ...
    ' $g = ', num2str(type_robust(2)),', w = ',num2str(type_robust(3)),'$']);


