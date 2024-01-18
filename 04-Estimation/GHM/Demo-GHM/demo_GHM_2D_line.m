%% Demo nonlinear GHM without constraints
%
% Wolfgang Förstner 6/6/2016
% wfoerstn@uni-bonn.de
%
% g = ^l2 - a ^l1 -b = 0 straight line  y_i - (a x_i + b) = 0
%
clearvars
close all

addpath(genpath('../../../General-Functions/'));
addpath(genpath('../Functions-GHM/'));

% initialization of random numbers
% = 0 CPU-time dependent
% > 0 fixed
%init_rand  = 0;
init_rand   = 15;     % standard = 15
init_rand = init_rand_seed(init_rand);

disp(' =============================================')
disp(' ----------  Demo GHM fitting line -----------')
disp(' ----------- y_i - (a x_i + b) = 0 -----------')
disp(' ---------------------------------------------')
disp(['Seed for random numbers                              = ' ,num2str(init_rand)]);

% choose type of covariance matrix
% cov_type = 0; % homogeneous, isotropic, uncorrelated
% cov_type = 1; % homogeneous, anisotropic, uncorrelated
cov_type = 2; % inhomogeneous, anisotropic, uncorrelated
%cov_type = 3; % random, correlated

% choose type of approximate values
xa = [0.8,1.0]';
% appr_type=0;
appr_type = 1;

%% generate data

% line parameters
a = 1;
b = 2;
disp(['true parameters  [a, b]                              = [', num2str(a),',',num2str(b),']'])

% number of observations
N = 4; .... cannot be changed
    
% true observations
l1_t  = randn(N,1);
l2_t  = a*l1_t + b;
% covariance matrix

switch cov_type
    case 0 % l1 and l2 independent, isotropic
        sigma_1 = 0.03;
        Cov_ll = [sigma_1^2*eye(N), zeros(N);zeros(N),sigma_1^2*eye(N)];
        % standard deviations
        disp(['Standard deviations for x_i and y_i                  =  ', num2str(sigma_1)])
        
    case 1
        % l1 and l2 independent, isotropic
        % standard deviations
        sigma_1 = 0.02;
        sigma_2 = 0.05;
        Cov_ll = [sigma_1^2*eye(N), zeros(N);zeros(N),sigma_2^2*eye(N)];
        disp(['Standard deviations for x_i and y_i                  =  ', num2str(sigma_1),',',num2str(sigma_2)])
        
    case 2
        % l1 and l2 independend, non isotropic
        var=[1,2,3,4,5,6,7,8];
        sigma_1=0.02;
        Cov_ll = diag(var)*sigma_1^2;
        disp(['Standard deviations for observations                 =  ', num2str(sqrt(var))])
    case 3
        % fully correlated
        sigma_1=0.01;
        A = randn(8);
        Cov_ll=(3*eye(8)+A*A')*sigma_1^2;
        disp('Standard deviations for observations                 =  random')
end

% all observations
lvt = [l1_t;l2_t];
% noisy observations
lv = rand_gauss(lvt,Cov_ll,1);
l1 = lv(1:N);
l2 = lv(N+1:2*N);


ScrS = plot_init;
figure('Color','w','Position',[100 100  ScrS+[ -600 -300]])
hold on
plot(l1,l2,'.r','MarkerSize',5)
for n=1:N
    ell.cov  = Cov_ll([n,n+4],[n,n+4]);
    ell.mean = [l1(n);l2(n)];
    plot_ellipse(ell,'-b');
end
title('Given points with standard ellipses, fitted line')
axis equal

r_U     = 1;                         % we are interested in the slope only
%% determine approximate values
if appr_type ==1 % estimate approximate values
    x0    = mean([l1,l2])';
    C0    = cov([l1,l2]);
    phi   = 0.5*atan2(2*C0(1,2),C0(2,2)-C0(1,1));
    xa(1) = tan(phi);
    xa(2) = x0(2) - x0(1) * xa(1); % xa(1) = tan(phi) = (x0(2)-xa(2))/(x0(1)-0)
end
disp(['approximate values                                   = [', num2str(xa(1)),',',num2str(xa(2)),']'])

%% estimate

sx      = [sigma_1,sigma_1]/sqrt(N); % approxmate standard deviations
Tx      = 0.02;                     % theshold for convergence
maxiter = 10;                        % maximum number of iterations

[est_x,Cov_xx,sigma_0q,R,vv,zv,riv,nabla_lv,muv,muv1,uv1q,uv2]...
    = GaussHelmertModel(lv,Cov_ll,@cg_2D_line,xa,sx,Tx,maxiter,r_U);

%% provide result
disp(' ')
disp(['estimated parameters  [a, b]   = [', num2str(est_x(1)),',',num2str(est_x(2)),']'])
disp(['estimated sigms_0              =  ', num2str(sqrt(sigma_0q))])
disp(['Covariance matrix              = [',num2str(Cov_xx(1,:)),']'])
disp(['                               = [',num2str(Cov_xx(2,:)),']'])
disp(['Standard deviations for [a, b] =  ', num2str(sqrt(Cov_xx(1,1))),', ',num2str(sqrt(Cov_xx(2,2)))])

xmin=min(l1)-0.2;
xmax=max(l1)+0.2;
ymin=min(l2);
ymax=max(l2);

plot([xmin,xmax],[est_x(1)*xmin+est_x(2),est_x(1)*xmax+est_x(2)],'-r','LineWidth',2);
xlim([xmin-0.2,xmax+0.2])
ylim([ymin-0.2,ymax+0.2])