%% Fig. 15.18 page 758
% Robust estimation for tackling steps and break points
%
% Wolfgang Förstner 2014-08-06
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de



close all
clearvars
 
addpath(genpath('../../General-Functions/'));
addpath('Functions')
 

%% plot settings
ss = plot_init;

%% set parameters

% params for generating the data

init_rand   = 6;                % just changes the noise, may be changed for other example
Nd          = 20;               % third of length of profile
N           = 3*Nd;             % total number of points
dens        = 1.0;              % percentage of observations
sigma_e     = 0.15;             % process noise
sigma_n     = 0.9;              % observation noise
Poutlier    = 0.00;             % P(outlier) - fraction of outliers
Max_outlier = 25;               % max(outlier)
factor = 1.0;                   % vertical scaling of true signal 
                                % (for changing signal to noise ratio)
type_outlier = 0;               % type of outlier (0=symm, 1=asymm)

% params for estimation
type_robust  = [0,...           %  [0 = L1, 1 = Kraus,g_factor,w_factor]
                2,...           %  g
                2,...           %  w
                1];             %  [0: not robust, 
                                %   1: only points,
                                %   2: only dem]
                        
print_type = 0;
plot_type  = 0;
                        
if type_outlier ==0
    Niter = 1;
    if type_robust(4) > 0
        Niter = 6;
    end
    factor_sigma = 1;
else
    Niter = 3;
    if type_robust(1) == 0
        factor_sigma= 8;
    else
        factor_sigma= 1;
    end
end


%% initialize random number generation by fixed seed
init_rand_seed(init_rand);

%% generate profile
[x,y,out_in,select,xs,ys] = ...
    generate_step_corner(Nd,sigma_e,sigma_n,Poutlier,Max_outlier,type_outlier,dens,factor);

%% reconstruct profile non robust
type_robust(4) = 0;   % reconstruction non robust
[xest_nonrob,~,~,weigths_nonrob] = estimate_profile_robust...
    (N,select,ys,sigma_e/factor_sigma,sigma_n,Niter,type_outlier,...
    type_robust,print_type,plot_type);

%% reconstruct profile robust
type_robust(4) = 1;     % reconstruction robust, only points
[xest_robust,~,v,weights] = estimate_profile_robust...
    (N,select,ys,sigma_e/factor_sigma,sigma_n,Niter,type_outlier,...
    type_robust,print_type,plot_type);

%% plot
figure('name','Fig. 16.18: Robust estimation for tackling steps and break points',...
    'color','w', 'Position',[0.1*ss(1),0.2*ss(2),0.5*ss(1),0.7*ss(2)]);

subplot(2,1,1);hold on
plot(1:N,x,'--r','LineWidth',2)
plot(1:N,xest_nonrob,'-k','LineWidth',2)
plot(1:N,ones(N,1),'-k')
plot(select,ys,'.b','MarkerSize',15)
plot(2:N-1,weigths_nonrob(2:N-1)+1,'.k','Markersize',15)
ylim([0,N/3+8]); xlim([0,N+1]); axis equal
xlabel('$x$'); ylabel('$z$');
title('Fig. 16.18a: Without robust estimation');


subplot(2,1,2);hold on
plot(1:N,x,'--r','LineWidth',2)
plot(1:N,xest_robust,'-k','LineWidth',2)
plot(1:N,ones(N,1),'-k')
plot(select,ys,'.b','MarkerSize',15)
plot(2:N-1,weights(end-(N-3):end)+1,'.k','Markersize',15)
xlim([0,N+1]); ylim([0,N/3+8]); axis equal
xlabel('$x$'); ylabel('$z$');
title({'Fig. 15.18b: With robust estimation.'...
    ['Iterations $1$ to $3$ use the $L_{12}$-norm,'...
    ' iterations $4$ to $6$ use the exponential weight function']});






