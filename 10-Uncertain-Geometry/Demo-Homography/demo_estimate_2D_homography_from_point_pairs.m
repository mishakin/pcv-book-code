%% demo_estimate_2D_homography_from_point_pairs
%
% generate 3x3 grid of points
% transform them with a homography
% add noise, possibly correlated between correspoinding points
%
% estimate homography (algebraical and ML-estimation)
%
% check estimation result (algebraic and ML-estimation)
%
% show scatter plots for transformed points
%      with theoretical covariance matrices
%
% see PCV Sect. 10.6.3, Fig. 10.28
% choose correlation in lines 36 ff.
%
% Wolfgang Förstner 2015
% wfoerstn@uni-bonn.de 


close all
clearvars
%clear all

addpath(genpath('../../General-Functions'))

global plot_option
global print_option
global print_option_estimation

disp('=========================================')
disp('      Estimate 2D homography             ')
disp('-----------------------------------------')
sugr_INIT
ss = plot_init;

%% control parameters

% random seed
init_rand      = 1;
%rinit_rand      = 0;

% number of cases
%N_samples     = 300;
%N_samples     = 100;
N_samples     = 50;
%N_samples     = 1;
disp(['Number of samples                      : ', num2str(N_samples)])

sigma_x       = 0.001;
disp(['Standard deviation of directions       : ', num2str(sigma_x)])
% correlation between points
%rho           = -0.98;
%rho           =     0;
rho           = +0.98;

disp(['Correlation of homologeous points      : ', num2str(rho)])

% Threshold for iteration
T       = 10^(-8);
maxiter = 10;

%significance level
S   = 0.999;
Slo = (1-S)/2;
Sup = 1-Slo;

%% random number intialization
init_rand_seed(init_rand);

%% Generate points in unit square
sigma_x       = 0.001;                  % std. dev.
sigma_y       = sigma_x;                  % std. dev.
N             = 9;                      % number of lines
H_tilde_e     = [1   -0.2  2.6;...
    0.1  1.2 -0.3; ...
    0.25 0.2    1];          % true homography matrix
H_tilde       = sugr_Homography_2D(H_tilde_e);   % true homography
factor_p      = 100;
b_random      = 0;                    % random points in [-1,+1]^2


print_option_estimation=1;
if N_samples >1
    print_option_estimation=0;
end
% do not plot samples
plot_option = 2;
print_option = 2;
if N_samples > 1
    plot_option = 0;
    print_option = 0;
end


% generate true points
PPt = sugr_generate_true_2D_point_pairs_homography(H_tilde.H,N,b_random);


est_H_samples_a    = zeros(N_samples,9);
est_H_samples_ml   = zeros(N_samples,9);
sigma_0s           = zeros(N_samples,1);
figure('Color','w','Position',[20 ss(2)/3 ss(1)/3 ss(2)/2])
hold on

start=cputime;
for i=1:N_samples
    if N_samples < 100
        fprintf(num2str(i)),fprintf(',')
        if mod(i,10)==0
            disp(' ')
        end
    else
        if mod(i,10)==0
            fprintf(num2str(i)),fprintf(',')
        end
        if mod(i,100)==0
            disp(' ')
        end
    end
    %disp(['sample: ',num2str(i)])
    % perturb point pairs
    PP = sugr_perturb_2D_point_pairs(PPt,sigma_x,sigma_y,rho,factor_p);
    
    %% Estimate point algebraically
    est_H_a = sugr_estimation_algebraic_Homography_2D_from_point_pairs(PP);
    
    %% Estimate point maximum likelihood
    [est_H_ml,sigma_0,R] = sugr_estimation_ml_Homography_2D_from_point_pairs(PP,est_H_a.H,T,maxiter);
    
    est_H_samples_a(i,:)    = est_H_a.H(:);
    est_H_samples_ml(i,:)   = est_H_ml.H(:);
    sigma_0s(i)             = sigma_0;
end
disp(['CPU time for ', num2str(N_samples),' samples               : ',num2str(cputime-start)]);


%% evaluate result



if N_samples > 2
    
    %% plot transferred points in 3x3 grid
    
    k=0;
    xp=zeros(N*N_samples,1);
    yp=zeros(N*N_samples,1);
    for i=1:N_samples
        for n=1:N
            %       given noisy points
            %[xn,yn,Cxy]   = sugr_select_Point_Pair_2D(PP,n);
            xn = PP.h(n,1:3)';
            yn = PP.h(n,4:6)';
            xne = xn(1:2)/xn(3);
            % true points as reference
            %[xt,yt,Cxy] = sugr_select_Point_Pair_2D(PPt,n);
            xt = PPt.h(n,1:3)';
            yt = PPt.h(n,4:6)';
            xte = xt(1:2)/xt(3);
            yte = yt(1:2)/yt(3);
            % true point transformed with estimated homography algebraic
            %y_a = sugr_transform_with_Homography_2D(est_H,xt);
            est_H = reshape(est_H_samples_a(i,:),3,3);
            y_a = est_H * xt;
            %ye_a = y_a.h(1:2)/y_a.h(3);
            ye_a = y_a(1:2)/y_a(3);
            %
            k=k+1;
            xp(k) = yte(1)+factor_p*(ye_a(1)-yte(1));
            yp(k) = yte(2)+factor_p*(ye_a(2)-yte(2));
            %plot(yte(1)+factor_p*(ye_a(1)-yte(1)),yte(2)+factor_p*(ye_a(2)-yte(2)),'.b');
        end
    end
    figure(1)
    hold on
    plot(xp,yp,'.b');
    
    
    k=0;
    xp=zeros(N*N_samples,1);
    yp=zeros(N*N_samples,1);
    for i=1:N_samples
        for n=1:N
            %       given noisy points
            %[xn,yn,Cxy]   = sugr_select_Point_Pair_2D(PP,n);
            xn = PP.h(n,1:3)';
            yn = PP.h(n,4:6)';
            xne = xn(1:2)/xn(3);
            % true points as reference
            %[xt,yt,Cxy] = sugr_select_Point_Pair_2D(PPt,n);
            xt = PPt.h(n,1:3)';
            yt = PPt.h(n,4:6)';
            xte = xt(1:2)/xt(3);
            yte = yt(1:2)/yt(3);
            % true point transformed with estimated homography algebraic
            %y_a = sugr_transform_with_Homography_2D(est_H,xt);
            est_H = reshape(est_H_samples_ml(i,:),3,3);
            y_a = est_H * xt;
            %ye_a = y_a.h(1:2)/y_a.h(3);
            ye_a = y_a(1:2)/y_a(3);
            %
            k=k+1;
            xp(k) = yte(1)+factor_p*(ye_a(1)-yte(1));
            yp(k) = yte(2)+factor_p*(ye_a(2)-yte(2));
            %plot(yte(1)+factor_p*(ye_a(1)-yte(1)),yte(2)+factor_p*(ye_a(2)-yte(2)),'.b');
        end
    end
    figure('Color','w','Position',[100+ss(1)/3  ss(2)/3 ss(1)/3 ss(2)/2])
    hold on
    plot(xp,yp,'.b');
    
    
    %
    est_H_mean_ml = mean(est_H_samples_ml);
    est_H_cov_ml  = cov(est_H_samples_ml);
    est_H_mean_a = mean(est_H_samples_a);
    est_H_cov_a  = cov(est_H_samples_a);
    
    est_H_ml_emp   = sugr_Homography_2D(reshape(est_H_mean_ml,3,3),est_H_cov_ml);
    est_H_a_emp   = sugr_Homography_2D(reshape(est_H_mean_a,3,3),est_H_cov_a);
    
    
    
    figure(1)
    %sugr_plot_Homography_2D(est_H_a_emp,'ok','-k',1,3*factor_p,'H');
    sugr_plot_Homography_2D(est_H_a,'ok','-k',1,3*factor_p,'H');
    title('ALG: samples with theoretical CovM')
    axis equal
    figure(2)
    sugr_plot_Homography_2D(est_H_ml,'ok','-k',1,3*factor_p,'H');
    title('ML: samples with theoretical CovM')
    axis equal
    
    R = 2*N-8;
    
    % for reducing the homography parameters (spectral normalization)
    Jhr = sugr_get_Jacobian_Jrh_Homography_2D(H_tilde.H);
    % Algebraic
    check_estimation_result(R,zeros(8,1),est_H_a.Crr,ones(N_samples,1),...
        est_H_samples_a*Jhr',S,' ALG')
    
    
    % Maximum likelihood
    check_estimation_result(R,zeros(8,1),est_H_ml.Crr,sigma_0s.^2,...
        est_H_samples_ml*Jhr',S,' ML');
    
end

% save(strcat('test_sugr_homography_estimation_a_ml_',num2str(rho*1000,'%4d')))