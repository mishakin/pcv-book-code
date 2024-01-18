%% Fig. 2.9 page 51 and Fig 16.11 page 742 
% simulate homgeneous/inhomogeneous 1D gaussian process
%
% Wolfgang Förstner 2015
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

close all
clearvars

addpath(genpath('../General-Functions'));

%% set your parameters ------------------------------------------
NoSample = 3;       % number of smaples for the plot, <5
N = 300;            % number of grid points

% parameters for homogeneous  gaussian process
d0 = 20;            % reference distance
sigma = 1;          % standard deviation

% parameters for inhomogeneous  gaussian process
d0l = 10;               % left reference distance
sigma_l = 0.5;          % left std
d0r = 40;               % right reference distance
sigma_r = 3;            % right standard deviation
dc = 30;                % smearing constant

%---------------------------------------------------------------

%% prepare visualization
% line spec
col=['k','b','r','c','m'];
lin=['-','-','-','-','-'];

% get current screensize, for proper positioning of figures and set default
% plot settings
ss = plot_init;
  
%% homogeneous GP

% initiate covariance matrix
Sigma=zeros(N);  
% generate 300x300 covariance matrix
for n = 1:N
    for m = 1:N
        d = (n-m)/d0;
        Sigma(n,m) = sigma^2*exp(-d^2/2);
    end
end

% generate three samples
y = rand_gauss(zeros(N,1),Sigma,NoSample);

% plot samples
figure('name','homogeneous GP','color','w','Position',[0.1*ss(1),0.2*ss(2),0.35*ss(1),0.60*ss(2)]); hold on
for s = 1:NoSample
    plot(1:N,y(:,s),strcat(lin(s),col(s)),'LineWidth',2);
end
plot(1:N, -3*sigma*ones(1,N),'--k')
plot(1:N, +3*sigma*ones(1,N),'--k')
xlim([-25,N+25])
ylim([-4.5,+4.5]*sigma)
xlabel('$t$');ylabel('$x$')
title('Fig. 2.9: Samples of homogeneous Gaussian Process')


