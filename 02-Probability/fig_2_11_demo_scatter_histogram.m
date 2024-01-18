%% Fig. 2.11 page 56
% example random number generation, scatterplot and histogram
%
% Wolfgang Förstner 2015
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

close all

addpath(genpath('../General-Functions/'));

% get current screensize, for proper positioning of figures and set default
% plot settings
ss = plot_init;

%% initialize random number generation by fixed seed 
%  (0 would use a random initialization)
init_rand_seed(2);

%% 1D random variable following a univariate Gaussian distribution

% parameters
N = 225;         % number of samples
mu_y = 2;        % mean
sigma_y = 0.25;   % standarddeviation

%% a) generate sample

% generate standard normal samples ~N(0,1)
x = randn(N,1);
% transform to desired Gaussian N(mu,sigma)
y = mu_y*ones(N,1)+sigma_y*x;

% plot samples
figure('name','Sample','color','w','Position',[0.1*ss(1),0.55*ss(2),0.35*ss(1),0.35*ss(2)]); hold on
for n=1:N
    plot([y(n),y(n)],[0,0.5],'-k','LineWidth',1)
end
plot([mu_y-4*sigma_y,mu_y+4*sigma_y],[0,0],'-k','LineWidth',1);
plot([0,0],[0,1],'-k','LineWidth',1);
xlim([mu_y-5*sigma_y,mu_y+5*sigma_y]); ylim([-0.1,1.1]);
xlabel('$y$');
% title(['Fig. 2.11a A sample of $N = 255$ normally distributed random variables $y\sim N(2,0.25)$'])
title(['Fig. 2.11a A sample of $N = ',num2str(N),'$ normally distributed random variables $y\sim N('...
    ,num2str(mu_y),',',num2str(sigma_y),')$'])

%% b) plot histogram (number of bins = sqrt(N))
figure('name','Histogram','color','w','Position',[0.55*ss(1),0.55*ss(2),0.35*ss(1),0.35*ss(2)]); hold on;

N_bin = floor(sqrt(N));       % number of bins
[NN,r]=hist(y,N_bin);         % calculate histogram
bar(r,NN)                     % plot histogram

% plot expected Gaussian / expected number of samples
range = abs(r(N_bin)-r(1))*N_bin/(N_bin-1);     % range of histogram
plot(r,N_bin*range*normpdf(r,mu_y,sigma_y),'-r','LineWidth',4);     
xlim([mu_y-5*sigma_y,mu_y+5*sigma_y]);ylim([0,max(NN)*1.1]);
xlabel('$y$');ylabel('$\sim p_y(y)$')
title(['Fig. 2.11b Histogram of the sample with $',num2str(N_bin),'$ bins, overlayed with its probability density.'])

%% 2d random vector with mean [0,0]

% parameters
sigma_x = 4.9;   % std in x
sigma_y = 3.2;   % std in y
rho = 0.7;       % correlation
N = 500;         % number of samples

% covariance matrix
C = [sigma_x^2 sigma_x*sigma_y*rho; sigma_x*sigma_y*rho sigma_y^2];
% generate sample
x = rand_gauss([0,0]',C,N)';

% plot samples with standard ellipse
figure('name','2D Sample','color','w','Position',[0.25*ss(1),0.1*ss(2),0.35*ss(1),0.35*ss(2)]); hold on
% standard ellipse 
plot_ellipse(struct('mean',[0;0],'cov',C),'-k');
% threefold standard ellipse
plot_ellipse(struct('mean',[0;0],'cov',9*C),'-g');  
% x-axis
plot([-4*sigma_x,+4*sigma_x],[0,0],'-k')
% y-axis
plot([0,0],[-4*sigma_y,+4*sigma_y],'-k')
% samples
plot(x(:,1),x(:,2),'.k')

xlim([-5*sigma_x,+5*sigma_x]);ylim([-5*sigma_y,+5*sigma_y])
axis equal
title({'Fig. 2.11c Sample of 2d random vector','standard ellipse (black), threefold standard ellipse (green)'})
