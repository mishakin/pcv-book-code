%% Fig. 2.10 page 54 
% Examples for autoregressive processes
%
% Calling this function yields four random profiles
% with the same parameters as in Fig. 2.10
%
% Repeated trials show typical samples
%
% Wolfgang Förstner 2015
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

close all
addpath(genpath('../General-Functions/'));

%% parameters
N = 100;                  % number of points
sigma_n = 0;              % observational noise
q = 0.999;                % decay of autoregressive processes
ymin = 100000;
ymax = -100000;

%% prepare figure
% get current screensize, for proper positioning of figures and set default
% plot settings
ss = plot_init;
figure('name','Examples for autoregressive processes','color','w','Position',[0.1*ss(1),0.2*ss(2),0.6*ss(1),0.7*ss(2)]); hold on

%% AR(1)     
p = 1; % order of process

% example a)
sigma_e = 1.0; % process noise
x11 = generate_observed_ARp(N,p,q,sigma_e,sigma_n);

% example b)
sigma_e = 0.2; % process noise
x12 = generate_observed_ARp(N,p,q,sigma_e,sigma_n);

%% AR(2)
p = 2;  % order of process

% example c)
sigma_e = 1.0; % process noise
x21 = generate_observed_ARp(N,p,q,sigma_e,sigma_n);

% example d)
sigma_e = 0.2; % process noise
x22 = generate_observed_ARp(N,p,q,sigma_e,sigma_n);


%% visualisation
% joint scale
ymin=min([ymin,min(x11),min(x12),min(x21),min(x22)]);
ymax=max([ymax,max(x11),max(x12),max(x21),max(x22)]);

subplot(2,2,1);
plot(1:N,x11,'.k','LineWidth',2)
xlim([1,N]);ylim([ymin,ymax]);xlabel('$n$');ylabel('$x$');
title('Fig. 2.10a AR(1), $\sigma_e = 1.0$')

subplot(2,2,2);
plot(1:N,x12,'.k','LineWidth',2);
xlim([1,N]);ylim([ymin,ymax]);xlabel('$n$');ylabel('$x$');
title('Fig. 2.10b AR(1), $\sigma_e = 0.2$')

subplot(2,2,3);
plot(1:N,x21,'.k','LineWidth',2)
xlim([1,N]);ylim([ymin,ymax]);xlabel('$n$');ylabel('$x$');
title('Fig. 2.10c AR(2), $\sigma_e = 1.0$')

subplot(2,2,4);
plot(1:N,x22,'.k','LineWidth',2)
xlim([1,N]);ylim([ymin,ymax]);xlabel('$n$');ylabel('$x$');
title('Fig. 2.10d AR(2), $\sigma_e = 0.2$')


