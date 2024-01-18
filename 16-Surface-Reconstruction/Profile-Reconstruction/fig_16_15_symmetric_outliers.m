%% Fig. 16.15 page 755
% test profile smoothing with symmetric outliers
%
% Wolfgang Förstner 2014-08-06
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

addpath(genpath('../../General-Functions/'));
addpath('Functions');

close all

%% plot settings
ss = plot_init;

disp('----- eliminate symmetric outliers -----')

%% set parameters

% params for generating the data
init_rand = 8;     % may be changed for other example
N = 200;           % number of points
dens = 0.8;        % percentage of observed points in (0,1]
sigma_e = 0.5;     % process noise
sigma_n = 0.5;     % observation noise
Poutlier = 0.4;    % fraction of outliers
Max_outlier = 25;  % max number of outliers
type_outlier = 0;  % symmetric  %type_outlier = 1; % asymmetric

% params for estimation
type_robust  = [0,2,2,2]; %  
    % 1. [0 = L1, 1 = Kraus,
    % 2. g_factor,
    % 3. w_factor
    % 4. 0,1,2,3 = 00,01,10,11 robust for points and dem

print_type = 0;
plot_type  = 0;

if type_outlier ==0
    Niter = 6;
    factor_sigma = 1;
else
    Niter = 3;
    factor_sigma= 8;
end


%% initialize random number generation by fixed seed
init_rand_seed(init_rand);

%% generate profile
[x,y,out_in,select,xs,ys] = ...
    generate_observed_AR2(N,sigma_e,sigma_n,Poutlier,Max_outlier,type_outlier,dens);

%% reconstruct profile
[xest,A,ver,weights,Cov] = ...
    estimate_profile_robust...
    (N,select,ys,sigma_e/factor_sigma,sigma_n,Niter,type_outlier,...
    type_robust,print_type,plot_type);


figure('name','Fig. 16.15 Reconstruction of profile with outliers','color','w',...
    'Position',[0.2*ss(1),0.2*ss(2),0.5*ss(1),0.5*ss(2)]);
hold on;
plot(1:N,x,'--r','LineWidth',3)
plot(1:N,xest,'-k','LineWidth',3)
plot(select,ys,'.b','MarkerSize',20)

title = (['Fig. 15.15: $N = ',num2str(N),'$, density $= ',num2str(dens),...
    '$, $s_e = ',num2str(sigma_e),'$, $s_n = ',num2str(sigma_e),...
    '$, $P($outlier$) = ',num2str(Poutlier),'$, $\max($outlier$) = ',num2str(Max_outlier),...
    '$, Niter$ = ',num2str(Niter),'$']);
xlabel('$x$');ylabel('$z$');

