%% Demo estimate 3D projection from 2D point pairs
%
% Model y = c(PX), perspective camera
%
% uncertain homogeneous coordinates for y,
% uncertain homogeneous coordinates for X
%
% simulation: 3D points X around origin, no points at infinity
%
% Wolfgang Förstner 6/2017
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 04/18
% wenzel@igg.uni-bonn.de

clc
close all
clearvars

addpath(genpath('../../General-Functions'))

global plot_option
global print_option
global print_option_estimation

sugr_INIT
ss = plot_init;

disp('============================================================')
disp(' Demo estimation of projection matrix, perspective model    ')
disp('         with uncertain image and scene points              ')
disp('             (generalizing PCV Alg. 17)                     ')
disp('------------------------------------------------------------')

%% control parameters
% for testing the implementation
%N_samples     = 1000;
N_samples     = 50;
%N_samples     = 1;
disp(strcat('Number of samples                   : ' , num2str(N_samples)));
% number of points
N              = 50;
%N               = 15;
%N             = 7;
disp(strcat('Number of points                    : ' , num2str(N)));

c = 500;
disp(strcat('Principal distance                  : ' , num2str(c)));

% std. dev.
sigma_X       =  0.0001;                  % std. dev. in m
sigma_y       = c/10000;                  % std. dev. in pixels
disp(strcat('Std. dev. of scene points           : ' , num2str(sigma_X)));
disp(strcat('Std. dev. of image points [pixel]   : ' , num2str(sigma_y)));

% random number intialization
 init_rand      = 12;                  % random seed
%init_rand      =  21378330;             % random seed

% Significance level for tests
S   = 0.999;
Slo = (1-S)/2;
Sup = 1- Slo;

% Threshold for iteration
T        = 10^(-10);
maxiter  = 20;
factor_p = 100;                  % for plot

%% random number intialization
init_rand_seed(init_rand);

%% Generate points in box
b_random      = 1;                        % random points in [-1,+1]^2

if b_random == 1
    dX            = 1;
    dY            = 1;
    dZ            = 0.5;
    Z0            = 3.0;
    Z_tilde   = [dX/2;dY/2;Z0];
else
    dX=0;dY=0;dZ=0;
    Z_tilde   = [0;0;3];
end

% set R and K
R_tilde       = calc_Rot_rod([1,2,3]);
K_tilde       = diag([c,c,1]);

% true projection
P_tilde         = sugr_Projection_3D_2D(K_tilde,R_tilde,Z_tilde);

disp('true_Projection');
disp(num2str(P_tilde.P));
disp(' ')

true_p          = P_tilde.P(:);

print_option_estimation = 1;
if N_samples > 1
    print_option_estimation = 0;
end

% generally do not plot samples
plot_optinon = 2;
print_option = 2;
% except for one sample
if N_samples > 1
    plot_option = 0;
    print_option = 0;
end

%% calculate samples

% generate true points, common to all samples
[Xt,yt] = sugr_generate_true_2D_point_pairs_Projection_3D_2D...
    (P_tilde.P,N,dX,dY,dZ,b_random);

X_coord  = Xt.e;
xs_coord = yt.e;

% set arrays for evaluation
est_P_samples_a_r  = zeros(N_samples,11);
est_P_samples_ml_r = zeros(N_samples,11);
est_s0q_samples_ml = zeros(N_samples,1);
est_bias_a_r       = zeros(N_samples,1);
est_bias_ml_r      = zeros(N_samples,1);
mean_rel_s         = zeros(N_samples,1);

start = cputime;
for i = 1:N_samples
    % monitor samples
    monitor_samples(i,N_samples);
    
    if i==1
        plot_option = 2;                                                   
    end
    % perturb point pairs
    [X,y] = sugr_perturb_3D_2D_point_pairs(Xt,yt,sigma_X,sigma_y,factor_p);
    set(gcf,'Position',[20,ss(2)/2-100,ss(1)/2,ss(2)/2]);
    
    X_coord  = X.e;
    xs_coord = y.e;
    plot_option = 0;
    
    
    %% Estimate point algebraically
    est_P_a = sugr_estimation_algebraic_Projection_3D_2D_from_point_pairs(X,y);
    P_est_algebraically = est_P_a.P;
    P_true = P_tilde.P;
    
    %% Estimate point maximum likelihood
    [est_P_ml,sigma_0,~] = ...
        sugr_estimation_ml_Projection_3D_2D_from_point_pairs...
        (X,y,P_tilde.P,T,maxiter);
    
    % store sigma_0^2
    est_s0q_samples_ml(i) = sigma_0^2;
    est_P_ml.P;
    P_tilde.P;
    
    eigv                = abs(eig(est_P_ml.Crr));
    condition_number    = max(eigv)/min(eigv);
    mean_rel_std        = sqrt(trace(inv(est_P_ml.Crr)*est_P_a.Crr)/11);    
    mean_rel_s(i)       = mean_rel_std;
    
    [Ke,Re,Ze] = calc_KRZ_from_P(est_P_ml.P,1);
    [K,R,Z]    = calc_KRZ_from_P(P_tilde.P,1);
    
    if print_option > 0
        disp(['K_Diff = ',inv(K)*Ke-eye(3)])                                
        disp(['R_Diff = ',inv(R)*Re-eye(3)])    
        disp(['Z_Diff = ',Z-Ze]) 
    end
    
    % collect estimated parameters
    est_P_samples_a_r(i,:)  = null(true_p')'*est_P_a.P(:);
    est_P_samples_ml_r(i,:) = null(true_p')'*est_P_ml.P(:);
    est_bias_a_r(i)  = ...
        est_P_samples_a_r(i,:)*inv(est_P_a.Crr)*est_P_samples_a_r(i,:)';
    est_bias_ml_r(i) = ...
        est_P_samples_ml_r(i,:)*inv(est_P_ml.Crr)*est_P_samples_ml_r(i,:)';
    
end
disp(['CPU time for ', num2str(N_samples),' samples              : ',num2str(cputime-start)]);

%% Check of estimation

if N_samples > 12
    % --------------------------------------------------------------------
    disp('Check estimation        ')
    disp('------------------------')
    
    R = 2*N-11;
    
    % Algebraic
    check_estimation_result(R,zeros(11,1),est_P_a.Crr,...
        ones(N_samples,1),est_P_samples_a_r,S,' ALG')
    
    % Maximum likelihood
    check_estimation_result(R,zeros(11,1),est_P_ml.Crr,...
        est_s0q_samples_ml,est_P_samples_ml_r,S,' ML');
    set(gcf,'Position',[ss(1)/3,20,ss(1)/3,ss(2)/3]);
    
    % show additional information
    figure('Color','w','Position',[ss(1)/2,ss(2)/2-100,ss(1)/2,ss(2)/2]);
    
    % relative accuracy: lambda(empALG/theorML)
    % algebraic estimation
    mean_p_a_r = mean(est_P_samples_a_r);
    covm_p_a_r = cov(est_P_samples_a_r);
    CovM_a_theor_r = est_P_a.Crr;
    
    % ml-estimation
    mean_p_ml_r = mean(est_P_samples_ml_r);
    covm_p_ml_r = cov(est_P_samples_ml_r);
    CovM_ml_theor_r = est_P_ml.Crr;
    %
    
    disp(['Mean rel. std.dev [sqrt(trace(C_a/C_ml)/11)] : ',...
        num2str(mean(mean_rel_s))]);
    
    % plots
    N_bin = ceil(sqrt(N_samples));
    
    subplot(2,2,1)
    hold on
    hist(est_s0q_samples_ml,N_bin)
    [a,r] = hist(est_s0q_samples_ml,N_bin);
    range = abs(r(N_bin)-r(1))*N_bin/(N_bin-1);     % range of histogram
    % expected number of samples
    plot(r,N_bin*range*fpdf(r,2*N-11,100000),'-r','LineWidth',4);
    title('estimated variance factor')
    
    subplot(2,2,2)
    hist(mean_rel_s,ceil(sqrt(N_samples)))
    title('mean relative std ALG/ML of P')
    
    subplot(2,2,3)
    hold on
    N_bin = ceil(sqrt(N_samples));
    [aa,r] = hist(est_bias_a_r,N_bin);
    hist(est_bias_a_r,N_bin)
    range = abs(r(N_bin)-r(1))*N_bin/(N_bin-1);     % range of histogram
    % expected number of samples
    plot(r,N_bin*range*chi2pdf(r,11),'-r','LineWidth',4);
    title('Mahalanobis d. of $P_{\mbox{ALG}}$')
    
    subplot(2,2,4)
    hold on
    N_bin = ceil(sqrt(N_samples));
    [aml,r] = hist(est_bias_ml_r,N_bin);
    hist(est_bias_ml_r,N_bin)
    range = abs(r(N_bin)-r(1))*N_bin/(N_bin-1);     % range of histogram
    % expected number of samples
    plot(r,N_bin*range*chi2pdf(r,11),'-r','LineWidth',4);
    title('Mahalanobis d. of $P_{\mbox{ML}}$')  

    
end

