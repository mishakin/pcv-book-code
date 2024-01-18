%% Fig. 16.14 page 754 
% compare different sampling distances
%
% Wolfgang Förstner 2014-08-06
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

addpath(genpath('../../General-Functions/'));
addpath('Functions')

close all

%% plot settings
ss = plot_init;

disp('..... compare different sampling distances -----')

%% set parameters

% params for generating the data
init_rand = 8;            % may be changed for other example
factorN = 80;             % ratio between grid distances
sample = [2,3,4,10,12];   % given positions in units of factorN
N = factorN*max(sample+1);% grid (corresponds to grid distance 1)
sample = sample*factorN;  % given positions in units 1
dens = 1.0;               % percentage of observed points in (0,1]
sigma_e = 0.5;            % process noise
sigma_n = 4;              % observation noise
factor_sigma = 1;         % factor for sigma_n
Poutlier = 0.0;           % fraction of outliers
Max_outlier = 25;         % max(outlier)
type_outlier = 0;         % symmetric   |  = 1; % asymmetric, following Kraus

% params for estimation
type_robust  = [0,2,2,0]; %  [0 = L1, 1 = Kraus,g_factor,w_factor]
print_type   = 0;
plot_type    = 0;

Niter = 0;

%% initialize random number generation by fixed seed
init_rand_seed(init_rand);

%% generate profile
[x,y,out_in] = ...
    generate_observed_AR2(N,sigma_e,sigma_n,Poutlier,...
    Max_outlier,type_outlier,dens);
 
% subsample
Nsample = length(sample);
select = zeros(Nsample,1);
xs = zeros(Nsample,1);
ys = zeros(Nsample,1);
for i=1:Nsample
    select(i) = sample(i);
    xs(i) = x(sample(i));
    ys(i) = y(sample(i));
end

%% reconstruct profile dense
[xest,A,ver,weights,Cov] = ...
    estimate_profile_robust...
    (N,select,ys,sigma_e,sigma_n*factor_sigma,Niter,type_outlier,type_robust,...
    plot_type,print_type);

%% reconstruct profile sparse
N_s = N/factorN;
select_s = select/factorN;
[xest_s,A_s,ver_s,weights_s,Cov_s] = ...
    estimate_profile_robust...
    (N_s,select_s,ys,sigma_e*factorN^(3/2),sigma_n*factor_sigma,Niter,...
    type_outlier,type_robust,plot_type,print_type);


%% plot
figure('name','Fig. 16.14: Test spacing','color','w',...
    'Position',[0.2*ss(1),0.2*ss(2),0.5*ss(1),0.5*ss(2)]);
hold on;
%plot(1:N,x,'--r','LineWidth',2)
plot(factorN-1:N,xest(factorN-1:N),'--k','LineWidth',2)
plot((1:N_s)*factorN,xest_s,'-b','LineWidth',2)
plot(select,ys,'.b','MarkerSize',15)

title(['Fig. 16.14: $N = ',num2str(N),'$ , high grid space $= ',num2str(factorN),...
    '$, $s_e = ',num2str(sigma_e),'$, $s_n = ',num2str(sigma_e),'$'],'FontSize',16);
xlabel('$x$');ylabel('$z$');

%% compare residuals
v1___v80___dx_div_sigma_n = ...
    [xest(sample)-ys,...
    xest_s(select_s)-ys,...
    (xest(sample)-xest_s(select_s))/sigma_n] %#ok<NOPTS>
