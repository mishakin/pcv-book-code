%% Fig. 16.9 page 740
% demo profile smoothing 
%
% Wolfgang Förstner 2014-10-07
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

addpath(genpath('../../General-Functions'))
addpath('Functions')

close all
clearvars

%% plot settings
ss = plot_init;

%% set parameters

% generation of points
N = 45;                % number of points

% estimation
Niter = 0;
type_robust = zeros(4,1);

% output
print_type = 0;
plot_type  = 0;
factor_sigma = 1.5;    % factor for displying std of estimated profile

%% ---------------------- LETS GO -----------------------------------------

disp('----- demo reconstruct profile with flattening and smoothing -----')

%% generate point
xs = [15,20,25,30]';
zs = [-1.8,-2,2,0]'*2.5;

%% reconstruct profile flat
sigma_e = 1.0;    % process noise
sigma_n = 0.25;   % observation noise

[zest_flat,tmp,tmp,tmp,Cov_flat] = estimate_profile_robust_flat...
    (N,xs,zs,sigma_e,sigma_n,Niter,0,type_robust,print_type);

% min_sigma = min(sqrt(diag(Cov_flat)))
% max_sigma = max(sqrt(diag(Cov_flat)))
% max_ratio = max(sqrt(diag(Cov_flat))/min_sigma)
% average_standard_deviation=sqrt(trace(Cov)/N)

%% reconstruct profile smooth
sigma_e = 1.0;     % process noise
sigma_n = 0.25;    % observation noise

[zest_smooth,tmp,tmp,tmp,Cov_smooth] = estimate_profile_robust...
    (N,xs,zs,sigma_e,sigma_n,Niter,0,type_robust,print_type,plot_type); %#ok<*ASGLU>

% min_sigma=min(sqrt(diag(Cov_smooth)))
% max_sigma=max(sqrt(diag(Cov_smooth)))
% max_ratio=max(sqrt(diag(Cov_smooth))/min_sigma)

%% plot
figure('name','Fig. 16.9 top: Profile smoothing','color','w',...
    'Position',[0.2*ss(1),0.2*ss(2),0.5*ss(1),0.5*ss(2)]);

% flat reconstruction
subplot(2,2,1); hold on
plot(1:N,zest_flat,'-k','LineWidth',2)
plot(xs,zs,'.b','MarkerSize',20)
axis equal
title('flat reconstruction','FontSize',16)
xlim([5,50]);ylim([-20,15])

% show precision
subplot(2,2,3); hold on
plot(xs,zs,'.b','MarkerSize',15)
plot(1:N,zest_flat,'-k','LineWidth',2)
plot(1:N,zest_flat+3*factor_sigma*sqrt(diag(Cov_flat)),'-r','LineWidth',1)
plot(1:N,zest_flat-3*factor_sigma*sqrt(diag(Cov_flat)),'-r','LineWidth',1)
title(['blow up factor for $\sigma_z = ',num2str(factor_sigma),'$'],'FontSize',16)
axis equal; xlim([5,50]); ylim([-20,15]);

% smooth reconstruction
subplot(2,2,2);hold on;
plot(1:N,zest_smooth,'-k','LineWidth',2)
plot(xs,zs,'.b','MarkerSize',20)
title('smooth reconstruction','FontSize',16)
axis equal;xlim([5,50]);ylim([-20,15]);

% show precision
subplot(2,2,4);hold on;
plot(xs,zs,'.b','MarkerSize',20)
plot(1:N,zest_smooth,'-k','LineWidth',2)
plot(1:N,zest_smooth+3*factor_sigma*sqrt(diag(Cov_smooth)),'-r','LineWidth',1)
plot(1:N,zest_smooth-3*factor_sigma*sqrt(diag(Cov_smooth)),'-r','LineWidth',1)
title(['blow up factor for $\sigma_z = ',num2str(factor_sigma),'$'],'FontSize',16)
axis equal; xlim([5,50]);ylim([-20,15]);



