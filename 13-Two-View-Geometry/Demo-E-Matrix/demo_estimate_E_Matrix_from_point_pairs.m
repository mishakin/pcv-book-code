%% test_estimate_E_Matrix_from_point_pairs
%
% direct solution following Sect. 13.3.2.3, p. 575 ff
%   with covariance matrix following Sect. 4.9.2.4
% iterative solution following Sect. 13.3.5, p. 585 ff
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 06/18
% wenzel@igg.uni-bonn.de

%clear all
clearvars
close all

disp(' ')
disp(' ')
disp('=================================================')
disp('---- Demo_estimate_E_Matrix_from_point_pairs ----')
disp('-------------------------------------------------')

global min_redundancy
global plot_option
global print_option
global print_option_estimation

addpath(genpath('../../General-Functions'))

sugr_INIT
ss = plot_init;

%% control parameters

%% random number intialization
init_rand      = 2;            
% init_rand      = 0;        
init_rand_seed(init_rand);

% number of corresponding points
% N             = 100;              
N             = 30;              
% N             = 15;                    

% number of samples for checking covariance matrix
%N_samples     = 1000;       
N_samples     = 50;         
%N_samples     = 1;       

%% standard deviations
% sigma_x       = 0.003;                  % std. dev. 
% sigma_y       = 0.003;                  % std. dev. 
%
% for Ladybug 3:  2*pi/(5*1400) rad/pixel, Lowe 0.3 pixel -> std radiant = 2.0944e-004
% sigma_x       = 0.00021;                  % std. dev. 
% sigma_y       = 0.00021;                  % std. dev. 
%
% for testing program
sigma_x       = 0.000001;                  % std. dev. 
sigma_y       = 0.000001;                  % std. dev. 

% correlation between points in left and right image
rho           = 0; %0.99;                     

% Threshold for iteration 
%(can be reduced to T = 0.1 and maxiter = 2 for practical purposes)
T       = 10^(-6);
maxiter = 5;
S = 0.999;

% in order to check also for low redundacies
min_redundancy = 2;

%% Generate points in unit cube

r_tilde       = [3.1,-2,-1.3];          %zeros(3,1); %
% r_tilde       = [0 0 0];                %zeros(3,1); %
R_tilde       = calc_Rot_r(r_tilde);    % true rotation matrix
b_tilde       = [0.2,1.0,1.35]';
% b_tilde       = [0 1 0]';
b_tilde       = b_tilde/norm(b_tilde);        % true basis
b_l_true      = 1;

% Representation of relative orientation without uncertainty
b_R = [b_tilde,R_tilde];
true_basis_Rotation_matrix = b_R                                           %#ok<*NOPTS>
% structure for relative orientation 
E_tilde       = sugr_E_Matrix(b_tilde,R_tilde);   % true basis, rotation
E_t           = E_tilde.E;
true_E_matrix = E_t

factor_p      = 12;                   % factor for plotting standard ellipses
b_random      = 1;                    % box for random points in [-1,+1]^3
Z0            = -1;                  % shift of bounding box (best to be negative)

print_option_estimation = 1;            % medium output
if N_samples >1 
    print_option_estimation = 0;
end

% do not plot samples
print_option = 2;
if N_samples > 1
    plot_option  = 0;
    print_option = 0;
end

%% generate true points
[PPt,Xt] = sugr_generate_true_2D_point_pairs_E_matrix(b_tilde*b_l_true,R_tilde,N,b_random,Z0);
set(gcf,'Name','Points and cameras','Position',[100 0.25*ss(2) 0.3*ss(1) 0.4*ss(2)])
set(gca,'CameraPosition',[-15, -8, 6])
%% initiate
est_E_samples_a  = zeros(N_samples,5);
est_bR_samples_a = zeros(N_samples,12);
est_E_samples_ml  = zeros(N_samples,5);
est_bR_samples_ml = zeros(N_samples,12);

total_time_a=0;
total_time_ml=0;

%% estimate
disp('Monitor samples:')
iii = 0;
s_a = zeros(1,N_samples);
s_ml = zeros(1,N_samples);
start = cputime;
for ii = 1:N_samples
    
    iii = iii+1;
    
    if N_samples < 100
        fprintf(num2str(ii)),fprintf(',')
        if mod(ii,10) == 0
            disp(' ')
        end
    else
        if mod(ii,10) == 0
            fprintf(num2str(ii)),fprintf(',')
        end
        if mod(ii,100) == 0
            disp(' ')
        end
    end
    
    %% perturb point pairs
    PP = sugr_perturb_2D_point_pairs_spherical(PPt,sigma_x,sigma_y,rho,factor_p);

    %% check
    % diag(PP.h(:,1:3)*calc_S(b_tilde)*R_tilde'*PP.h(:,4:6)')';
    
    %% Estimate E_matrix algebraically
    start = cputime;

    [est_E_a,sigma_0_a,error] = sugr_estimation_algebraic_E_Matrix_from_point_pairs(PP);
    if error > 0
        disp(['sample ',num2str(ii),' failed'])
        iii = iii-1;
    else
        check_algebraic = diag(PP.h(:,1:3)*est_E_a.E*PP.h(:,4:6)')';
        
        b_a_r                  = null(b_tilde')' * est_E_a.bR(:,1);
        [r_a,a_a]              = calc_r_phi_from_R(est_E_a.bR(:,2:4)*R_tilde');
        est_E_samples_a(ii,:)  = [b_a_r',r_a'*a_a];
        est_bR_samples_a(ii,:) = est_E_a.bR(:)';

        s_a(ii)  = 1;

        total_time_a = total_time_a + cputime-start;
        %% Estimate E_matrix maximum likelihood
        start = cputime;

        [est_E_ml,sigma_0_ml,R] = sugr_estimation_ml_E_Matrix_from_point_pairs(PP,est_E_a.bR,T,maxiter);
        total_time_ml           = total_time_ml + cputime-start;
        

        b_ml_r                  = null(b_tilde')' * est_E_ml.bR(:,1);
        [r_ml,a_ml]             = calc_r_phi_from_R(est_E_ml.bR(:,2:4)*R_tilde');
        est_E_samples_ml(ii,:)  = [b_ml_r',r_ml'*a_ml];
        est_bR_samples_ml(ii,:) = est_E_ml.bR(:)';

        s_ml(ii) = sigma_0_ml;
    end
end
disp(['CPU time for ', num2str(N_samples),' samples              : ',num2str(cputime-start)]);
disp(' ')

%%
rho_degree = 180/pi;
disp(strcat('sigma basis direction [mrad] = ',...
    num2str(1000*sqrt(diag(est_E_ml.CbRbR(1:2,1:2)))')));
disp(strcat('sigma rotation        [mrad] = ',...
    num2str(1000*sqrt(diag(est_E_ml.CbRbR(3:5,3:5)))')));
disp(strcat('sigma basis direction [degree] = ',...
    num2str(rho_degree*sqrt(diag(est_E_ml.CbRbR(1:2,1:2)))')));
disp(strcat('sigma rotation        [degree] = ',...
    num2str(rho_degree*sqrt(diag(est_E_ml.CbRbR(3:5,3:5)))')));

%% check algorithm
if N_samples > 5 
    
    % check of ALG 
    check_estimation_result(N-5,zeros(5,1),est_E_a.CbRbR,s_a.^2,...
        est_E_samples_a,S,' ALG')
   
    % ML 
    check_estimation_result(N-5,zeros(5,1),est_E_ml.CbRbR,s_ml.^2,...
        est_E_samples_ml,S,' ML');
    
    set(gcf,'Name','Histogram var.factors','Position',[0.4*ss(1) 0.25*ss(2) ss(1)/2 0.4*ss(2)])
end

total_CPUtime_alg_mle_seconds  = [total_time_a,total_time_ml]

CovM_dir = sugr_CovM_E_matrix(Xt, PP.Crr, b_tilde, R_tilde);
CovM_comparison_generalized_eigenvalues_should_be_one = ...
    eigs(inv(est_E_ml.CbRbR)*CovM_dir)'                                    %#ok<MINV>

disp('          ========= End of Demo ===========')
disp(' ')
disp(' ')
