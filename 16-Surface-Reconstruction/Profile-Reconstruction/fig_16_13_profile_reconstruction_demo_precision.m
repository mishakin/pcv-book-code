%% Fig 16.13 page 750
% profile recondstruction with precision
%
% Wolfgang Förstner 2014-08-06
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

function fig_16_13_profile_reconstruction_demo_precision()

addpath(genpath('../../General-Functions/'));
addpath('Functions')

close all

%% initialize random number generation by fixed seed
init_rand_seed(23);

%% plot settings
ss = plot_init;

%% set parameters
N = 200;                    % number of grid points
sigma_e=0.5;                % process noise
sigma_n=0.5;                % observation noise
factor_sigma = 1;           % factor for sigma_n
type_outlier = 0;           % symmetric
type_robust  = [0,0,0,0];   % not robust
Niter = 0;                  % not robust
print_type = 0;
plot_type  = 0;

%% generate profile
[x,y,select,xs,ys] = ...
    generate_observed_AR2_demo(N,sigma_e,sigma_n);

%% reconstruct profile
[xest,aa,ab,ac,Cov] = estimate_profile_robust...
    (N,select,ys,sigma_e/factor_sigma,sigma_n,Niter,...
    type_outlier,type_robust,...
    plot_type,print_type);

%% show precision

factor = 1;               % blow-up factor for error band

figure('name','Fig. 16.13: Profile Reconstruction and precision','color','w',...
    'Position',[0.2*ss(1),0.2*ss(2),0.5*ss(1),0.5*ss(2)]);
hold on;
% for n=1:N-1
%     plotrect([n,x(n)-3*factor*sqrt((Cov(n,n)))],...
%              [n+1,x(n+1)+3*factor*sqrt((Cov(n,n)))],'g' );
% end
xbandtop = x+3*factor*sqrt(diag(Cov));
xbandbottom = x-3*factor*sqrt(diag(Cov));
fill([(1:N)';(N:-1:1)'],[xbandtop;xbandbottom(end:-1:1)],[0.3,0.97,0.02],'FaceAlpha',0.5,'EdgeAlpha',0);
plot(1:N,xest,'-k','LineWidth',2)
plot(1:N,x,'--r','LineWidth',2)
plot(1:N,xbandtop,'-r','LineWidth',1)
plot(1:N,xbandbottom,'-r','LineWidth',1)
plot(select,ys,'.b','MarkerSize',15)
plot(1:N,xest,'.k','LineWidth',2)
title(['Fig 16.13: $N = ',num2str(N),'$, $s_e = ',num2str(sigma_e),'$, $s_n =',num2str(sigma_e),'$'])


return


%% generate_observed_AR2
% 
% N         = number of points
% sigma_e   = process noise
% sigma_n   = observation noise
%
% x,y       = true signal
% select    = indices for observed points
% xy,ys     = observed points

function [x,y,select,xs,ys] = ...
    generate_observed_AR2_demo(N,sigma_e,sigma_n)

x=zeros(N,1);
y=zeros(N,1);

for n = 3:N
    x(n) = 1.9998*x(n-1)-0.9999*x(n-2)+randn(1)*sigma_e;
    y(n) = x(n)+randn(1)*sigma_n;
end
for i=1:N
    y(i) = y(i) - (i-1)*(x(N)-x(1))/(N-1);
    x(i) = x(i) - (i-1)*(x(N)-x(1))/(N-1);
end

M = [6,12,24,81,95,124,138,176];
select = zeros(8,1);
xs = zeros(8,1);
ys = zeros(8,1);
for m = 1:8
    select(m) = M(m);
    xs(m) = x(M(m));
    ys(m) = y(M(m));
end

return
