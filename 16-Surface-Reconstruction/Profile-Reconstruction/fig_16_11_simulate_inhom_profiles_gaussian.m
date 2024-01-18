%% Fig. 16.11 on p. 742 
% simulate homgeneous/inhomogeneous 1D gaussian process
%
% Wolfgang Förstner 2015
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

close all
clearvars

addpath(genpath('../../General-Functions'));

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
  

%% inhomogeneous GP

% initiate covariance matrix
Sigma=zeros(N);  
% generate left upper covariance matrix
for n=1:N/2
    for m=1:N/2
        d=(n-m)/d0l;
        % smear sigma using sigmoid function
        sig = (sigma_r+(sigma_l-sigma_r)./(1+exp(-(n-N/2)./N.*dc)))*...
              (sigma_r+(sigma_l-sigma_r)./(1+exp(-(m-N/2)./N.*dc)));
        Sigma(n,m)=sig*exp(-d^2/2);
    end
end
% generate right lower covariance matrix
for n=N/2+1:N
    for m=N/2+1:N
        d=(n-m)/d0r;
        % smear sigma using sigmoid function
        sig = (sigma_r+(sigma_l-sigma_r)./(1+exp(-(n-N/2)./N.*dc)))*...
              (sigma_r+(sigma_l-sigma_r)./(1+exp(-(m-N/2)./N.*dc)));
        Sigma(n,m)=sig*exp(-d^2/2);
    end
end
% generate  right upper covariance matrix
for n=1:N/2
    for m=N/2+1:N
        % take mean distance
        d=abs(N/2+1/2-n)/d0l + abs(m-N/2-1/2)/d0r;
        % smear sigma using sigmoid function
        sig = (sigma_r+(sigma_l-sigma_r)./(1+exp(-(n-N/2)./N.*dc)))*...
              (sigma_r+(sigma_l-sigma_r)./(1+exp(-(m-N/2)./N.*dc)));
        Sigma(n,m) = sig*exp(-d^2/2);
        % fill symmetric part
        Sigma(m,n) = Sigma(n,m);
    end
end

% generate samples
y = rand_gauss(zeros(N,1),Sigma,NoSample);

% prepare plot
rang = N:-1:1;
figure('name','inhomogeneous GP','color','w','Position',[0.55*ss(1),0.2*ss(2),0.35*ss(1),0.60*ss(2)]);hold on
% plot samples: y + 3*sigma(smeared)
for s = 1:NoSample
    plot(1:N,y(rang,s)'+...
        3*(sigma_r+...
        (sigma_l-sigma_r)./(1+exp(-(rang-N/2)./N*dc))),...
        strcat(lin(s),col(s)),'LineWidth',3);
end
xlabel('$x$');ylabel('$z$')
title('Fig. 16.11: Samples of inhomogeneous profiles')
