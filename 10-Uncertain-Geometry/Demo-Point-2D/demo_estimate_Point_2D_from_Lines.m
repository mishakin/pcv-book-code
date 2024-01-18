%% demo_estimate_Point_2D_from_Lines
%
% see Sect. 10.6.2.1, Figs. 10.24-10.26
%
% four types of estimation of intersection point (vanishing point)
%       ALGe algebraic Euclideanly
%       ALGs algebraic spherically
%       SSDe geometric Euclideanly
%       MLEe ML Euclideanly
%       MLEs ML spherically
%
% comparison of empirical and theoretical covariance matrices
% 
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 06/18
% wenzel@igg.uni-bonn.de

addpath(genpath('../../General-Functions'))

global min_redundancy
global plot_option
global print_option
global print_option_estimation

close all
% clear all
clearvars
clc

sugr_INIT
ss = plot_init;
min_redundancy=100000;

disp('=====================================================================')
disp('   Point estimation (geometric, ALG, ML; Euclideanly, spherically)   ')
disp('                    Figures 10.24 - 10.26                            ')
disp('---------------------------------------------------------------------')

%% Control parameters
N_samples     = 100;                  % number of cases (for first trial)
%N_samples     = 1000;                  % number of cases
disp(strcat('Number of samples            : ' , num2str(N_samples)));

% sample sizes
N             = 50;                   % number of lines per case
disp(strcat('Number of lines              : ' , num2str(N)));

% parameters for line generation
d             = 1.0;                  % size of square [-d,+d]^2 image wrt c=1
resolution    = 1000;                 % number of pixels per d


disp(strcat('Size of image                : ' , num2str(resolution)));

%
lengthq       = 200;                  % average length of line segments
disp(strcat('Average length of lines      : ' , num2str(lengthq)));
sigma         = 1.0;                  % std. dev. of pixels on line
disp(strcat('Std. dev. of pixels on lines : ' , num2str(sigma)));

% initial seed
init_rand      = 11;                  % random seed

% type of simulation
sim_type = 0; % close unknown point  -->  Figs. 10.24-10.25
%sim_type = 1; % far unknown point     -->  Fig.  10.26

% factor for plotting the standard ellipses
factor_p      = 1000;   

% rigorous variance propagation for observed line segment
%rigorous      = 0;
rigorous      = 1;

%% initialize random seed
init_rand_seed(init_rand);

%% Generate lines           

if N_samples >1 
    print_option_estimation=0;
end
%
switch sim_type
    case 0
        x_tilde_e      =[1.3;1.7];
    case 1
        x_tilde_e     = [0,200]';
end
disp(strcat('True point                   : [',num2str(x_tilde_e(1)),',',num2str(x_tilde_e(2)),']'));

%
x_tilde       = sugr_Point_2D(x_tilde_e,10^-10*eye(2));   % true point

% -------------------------------------------------------------------------
% generate true lines
lines = sugr_generate_true_Lines_2D_one_Point(x_tilde.h,d,N,lengthq,resolution);

% show representatative result
plot_option = 2;                                                           %#ok<NASGU>
print_option = 2;                                                          %#ok<NASGU>

if sim_type == 0    
    
    figure('name','estimate intersection points','Color','w','Position',[10,50,0.5*ss(1),0.85*ss(2)])
    subplot(3,2,1)
    hold on
    sugr_perturb_Lines_2D(x_tilde.h,d,lines,resolution,sigma,rigorous);    
    title('Original line segments')
    
    subplot(3,2,2)
    hold on
    sugr_perturb_Lines_2D(x_tilde.h,d,lines,resolution,sigma,rigorous);
    
    subplot(3,2,3)
    hold on
    sugr_perturb_Lines_2D(x_tilde.h,d,lines,resolution,sigma,rigorous);
        
    subplot(3,2,4)
    hold on
    sugr_perturb_Lines_2D(x_tilde.h,d,lines,resolution,sigma,rigorous);
        
    subplot(3,2,5)
    hold on
    sugr_perturb_Lines_2D(x_tilde.h,d,lines,resolution,sigma,rigorous);
    
    subplot(3,2,6)
    hold on
    sugr_perturb_Lines_2D(x_tilde.h,d,lines,resolution,sigma,rigorous);
    
end

% do not plot samples in the following
plot_option  = 0;
print_option = 0;
% -------------------------------------------------------
est_x_samples_aE = zeros(N_samples,2);
sigma_2_aE_D     = zeros(N_samples,1);
sigma_2_aE_O     = zeros(N_samples,1);

est_x_samples_a = zeros(N_samples,2);
sigma_2_a_D     = zeros(N_samples,1);
sigma_2_a_O     = zeros(N_samples,1);

est_x_samples_g = zeros(N_samples,2);
sigma_2_g_D     = zeros(N_samples,1);
sigma_2_g_O     = zeros(N_samples,1);

est_x_samples_g_x = zeros(N_samples,2);
sigma_2_g         = zeros(N_samples,1);

est_x_samples_ml = zeros(N_samples,2);
sigma_2_ml        = zeros(N_samples,1);
sigma_2_ml_D     = zeros(N_samples,1);
sigma_2_ml_O     = zeros(N_samples,1);

est_x_samples_mH = zeros(N_samples,2);
sigma_2_mH        = zeros(N_samples,1);
sigma_2_mH_D     = zeros(N_samples,1);
sigma_2_mH_O     = zeros(N_samples,1);

start = cputime;
for n = 1:N_samples
     if N_samples < 100
        fprintf(num2str(n)),fprintf(',')
        if mod(n,10) == 0
            disp(' ')
        end
    else
        if mod(n,10) == 0
            fprintf(num2str(n)),fprintf(',')
        end
        if mod(n,100) == 0
            disp(' ')
        end
    end
    % perturb lines
    lstruct = sugr_perturb_Lines_2D(x_tilde.h,d,lines,resolution,sigma,rigorous);
    
    
    %% Estimate point algebraically - Euclideanly
    [est_p_aE,sigma_0_est_DE,sigma_0_est_OE] = ...
        sugr_estimation_algebraic_Point_2D_from_Lines(lstruct,0);
    est_x_aE_e = est_p_aE.h(1:2)'/est_p_aE.h(3);
    est_x_samples_aE(n,:) = est_x_aE_e;
    sigma_2_aE_D(n) = sigma_0_est_DE^2;
    sigma_2_aE_O(n) = sigma_0_est_OE^2;    
    
    if sim_type == 0
        plot([-2,-2,2,2,-2],[-2,2,2,-2,-2],'-b')
        
        subplot(3,2,2)
        hold on
        plot(x_tilde_e(1)+factor_p*(est_x_aE_e(1)-x_tilde_e(1))/3,...
            x_tilde_e(2)+factor_p*(est_x_aE_e(2)-x_tilde_e(2))/3,'.g');
    end
    
    %% Estimate point algebraically - spherically
    [est_p_a,sigma_0_est_D,sigma_0_est_O] = ...
        sugr_estimation_algebraic_Point_2D_from_Lines(lstruct,1);
    est_x_a_e = est_p_a.h(1:2)'/est_p_a.h(3);
    est_x_samples_a(n,:) = est_x_a_e;
    sigma_2_a_D(n) = sigma_0_est_D^2;
    sigma_2_a_O(n) = sigma_0_est_O^2;    
    
    if sim_type == 0
        plot([-2,-2,2,2,-2],[-2,2,2,-2,-2],'-b')        
        
        subplot(3,2,3)
        hold on
        plot(x_tilde_e(1)+factor_p*(est_x_a_e(1)-x_tilde_e(1))/3,...
            x_tilde_e(2)+factor_p*(est_x_a_e(2)-x_tilde_e(2))/3,'.g');
    end
    
    %% Estimate point geometrically
    [est_p_g,sigma_0_g,xe] = sugr_estimation_geometric_Point_2D_from_Lines(lstruct);
    est_x_g_e              = est_p_g.h(1:2)'/est_p_g.h(3);
    est_x_samples_g(n,:)   = est_x_g_e;
    est_x_samples_g_x(n,:) = xe;
    sigma_2_g(n)           = sigma_0_g^2;
    
    if sim_type == 0
        plot([-2,-2,2,2,-2],[-2,2,2,-2,-2],'-b')        
        
        subplot(3,2,4)
        hold on
        plot(x_tilde_e(1)+factor_p*(est_x_g_e(1)-x_tilde_e(1))/3,...
            x_tilde_e(2)+factor_p*(est_x_g_e(2)-x_tilde_e(2))/3,'.g');
    end
    
    %% estimate ML
    T = 0.01;
    maxiter = 3;
    [est_p_ml,sigma_0_ml,~] = sugr_estimation_ml_Point_2D_from_Lines(lstruct,est_p_a.h,T,maxiter);
    [est_x_ml_e,Cxx_ml]     = sugr_get_Euclidean_Point_2D(est_p_ml);
    est_x_samples_ml(n,:)   = est_x_ml_e;
    sigma_2_ml(n)           = sigma_0_ml^2;
    
    if sim_type == 0
        plot([-2,-2,2,2,-2],[-2,2,2,-2,-2],'-b')
        
        subplot(3,2,5)
        hold on
        plot(x_tilde_e(1)+factor_p*(est_x_ml_e(1)-x_tilde_e(1))/3,...
            x_tilde_e(2)+factor_p*(est_x_ml_e(2)-x_tilde_e(2))/3,'.g');
    end
    %% estimate ML-Hessian
    T = 0.01;
    maxiter = 3;
    [est_p_mH,sigma_0_mH,R] = sugr_estimation_ml_Point_2D_from_Lines_Hessian(lstruct,est_p_a.h,T,maxiter);
    [est_x_mH_e,Cxx_mH]  = sugr_get_Euclidean_Point_2D(est_p_mH);
    est_x_samples_mH(n,:) = est_x_mH_e;
    sigma_2_mH(n) = sigma_0_mH^2;    
    
    if sim_type == 0
        plot([-2,-2,2,2,-2],[-2,2,2,-2,-2],'-b')        
        
        subplot(3,2,6)
        hold on
        plot(x_tilde_e(1)+factor_p*(est_x_mH_e(1)-x_tilde_e(1))/3,...
             x_tilde_e(2)+factor_p*(est_x_mH_e(2)-x_tilde_e(2))/3,'.g');        
        
        plot([-2,-2,2,2,-2],[-2,2,2,-2,-2],'-b')
    end
end
disp(['CPU time for ', num2str(N_samples),' samples              : ',num2str(cputime-start)]);

%% Analyse samples
[x_aE_e,Cxx_aE_e] = sugr_get_Euclidean_Point_2D(est_p_aE);
[x_a_e, Cxx_a_e]  = sugr_get_Euclidean_Point_2D(est_p_a);
[x_g_e, Cxx_g_e]  = sugr_get_Euclidean_Point_2D(est_p_g);
[x_ml_e,Cxx_ml_e] = sugr_get_Euclidean_Point_2D(est_p_ml);
[x_mH_e,Cxx_mH_e] = sugr_get_Euclidean_Point_2D(est_p_mH);

if N_samples > 2
    %% determine empirtical means and covariance matrices
    est_x_mean_aE  = mean(est_x_samples_aE);
    est_x_cov_aE   = cov(est_x_samples_aE);
    est_x_mean_a   = mean(est_x_samples_a);
    est_x_cov_a    = cov(est_x_samples_a);
    est_x_mean_g   = mean(est_x_samples_g);
    est_x_cov_g    = cov(est_x_samples_g);
    est_x_mean_ml  = mean(est_x_samples_ml);
    est_x_cov_ml   = cov(est_x_samples_ml);
    est_x_mean_mH  = mean(est_x_samples_mH);
    est_x_cov_mH   = cov(est_x_samples_mH);
    
    mean_sigma_2_a_OE = mean(sigma_2_aE_O);
    mean_sigma_2_a_DE = mean(sigma_2_aE_D);
    mean_sigma_2_a_O  = mean(sigma_2_a_O);
    mean_sigma_2_a_D  = mean(sigma_2_a_D);
    mean_sigma_2_g    = mean(sigma_2_g);
    mean_sigma_2_ml_Euclidean   = mean(sigma_2_ml);
    mean_sigma_2_ml_spherical   = mean(sigma_2_mH);
    
    %%
    est_p_aE_emp    = sugr_Point_2D(est_x_mean_aE', est_x_cov_aE);
    est_p_a_emp     = sugr_Point_2D(est_x_mean_a', est_x_cov_a);
    est_p_g_emp     = sugr_Point_2D(est_x_mean_g', est_x_cov_g);
    est_p_ml_emp    = sugr_Point_2D(est_x_mean_ml',est_x_cov_ml);
    est_p_mH_emp    = sugr_Point_2D(est_x_mean_mH',est_x_cov_mH);    
    
    if sim_type == 0
        subplot(3,2,1)
        if sim_type == 0
            xlim([-2.2,+2.5])
            ylim([-2.2,2.9])
        end
        
        subplot(3,2,2)
        
        title('ALGe emp (green) - ALGe theor (blue)')
        hold on
        sugr_plot_Point_2D(est_p_aE_emp,'ok','-g',2,factor_p);
        hold on
        sugr_plot_Point_2D(est_p_aE,'ok','-b',1,factor_p);
        axis equal
        if sim_type == 0
            xlim([-2.2,+2.5])
            ylim([-2.2,2.9])
        end
        
        subplot(3,2,3)
        
        title('ALGs')
        hold on
        %sugr_plot_Point_2D(est_p_a_emp,'ok','-g',2,factor_p);
        hold on
        sugr_plot_Point_2D(est_p_a,'ok','-k',1,factor_p);
        axis equal
        if sim_type == 0
            xlim([-2.2,+2.5])
            ylim([-2.2,2.9])
        end
        
        subplot(3,2,4)
        title('SSDe')
        hold on
        %sugr_plot_Point_2D(est_p_g_emp,'ok','-g',2,factor_p);
        hold on
        %sugr_plot_Point_2D(est_p_g,'ok','-b',1,factor_p);
        axis equal
        if sim_type == 0
            xlim([-2.2,+2.5])
            ylim([-2.2,2.9])
        end
        
        subplot(3,2,5)
        title('MLEe')
        % sugr_plot_Point_2D(est_p_mH_emp,'ok','-g',2,factor_p);
        hold on
        sugr_plot_Point_2D(est_p_mH,'ok','-k',1,factor_p);
        axis equal
        if sim_type == 0
            xlim([-2.2,+2.5])
            ylim([-2.2,2.9])
        end        
        
        subplot(3,2,6)
        title('MLEs (green) - ALGs (dashed)')
        sugr_plot_Point_2D(est_p_a,'ok','--k',1,factor_p);
        hold on
        sugr_plot_Point_2D(est_p_ml,'ok','-k',1,factor_p);
        axis equal
        if sim_type == 0
            xlim([-2.2,+2.5])
            ylim([-2.2,2.9])
        end
    end
    
% % distance between the covariance matrices
% distance_a_emp= exp(sqrt(sum(log(eigs(inv(Cxx_a_e)*est_x_cov_a)).^2)/2)/2)
% distance_g_emp= exp(sqrt(sum(log(eigs(inv(Cxx_g_e)*est_x_cov_g)).^2)/2)/2)
% 
% distance_a_ml = exp(sqrt(sum(log(eigs(inv(Cxx_ml_e)* Cxx_a_e )).^2)/2)/2)
% distance_g_ml = exp(sqrt(sum(log(eigs(inv(Cxx_ml_e)* Cxx_g_e )).^2)/2)/2)
end

%% final distribution of points

% algebraic solution
if N_samples > 2
%% show histograms
    figure('name','Histogramms','Color','w','Position',[ss(1)/2,0.3*ss(2),0.5*ss(1),0.6*ss(2)])   
    N_bin = ceil(sqrt(N_samples));
    
    subplot(3,2,2)
    hist(sigma_2_a_D,N_bin);
    title('Variance ALGe')
    
    subplot(3,2,3)
    hist(sigma_2_aE_D,N_bin);
    title('Variance ALGs')
    
    subplot(3,2,4)
    hist(sigma_2_g,N_bin);
    title('Variance SSDe')
    
    subplot(3,2,5)
    hold on
    hist(sigma_2_mH,N_bin);
    [~,r] = hist(sigma_2_mH,N_bin);
    %r=2*(1:ceil(sqrt(N_samples))-1/2)/N_bin;
    range = abs(r(N_bin)-r(1));
    plot(r,N_bin*range*fpdf(r,R,10000),'-r','LineWidth',2);
    
    title('Variance factor MLEe')
    
    subplot(3,2,6)
    hold on
    hist(sigma_2_ml,N_bin);
    [NN,r] = hist(sigma_2_ml,N_bin);
    %r=2*(1:ceil(sqrt(N_samples))-1/2)/N_bin;
    range = abs(r(N_bin)-r(1));
    plot(r,N_bin*range*fpdf(r,R,10000),'-r','LineWidth',2);
    
    title('Variance factor MLEs')

    %%  show 
    figure('name','Distribution of points','Color','w','Position',[100,20,0.6*ss(1),0.82*ss(2)])
    plot_option=1;
    subplot(3,2,1)
    hold on
    title('Straight line segments')
    sugr_perturb_Lines_2D(x_tilde.h,d,lines,resolution,sigma,rigorous);

    subplot(3,2,2)
    hold on
    title('ALGe')
    for n=1:N_samples
        plot(est_x_mean_aE(1)+(est_x_samples_aE(n,1)-est_x_mean_aE(1)),...
            est_x_mean_aE(2)+(est_x_samples_aE(n,2)-est_x_mean_aE(2)),'.k');
    end
    sugr_plot_Point_2D(est_p_aE_emp,'ok','--k',1,3);

    plot(est_x_mean_aE(1),est_x_mean_aE(2),'ow','MarkerSize',12,'LineWidth',7);
    plot(est_x_mean_aE(1),est_x_mean_aE(2),'ok','MarkerSize',12,'LineWidth',3);

    plot(x_tilde.h(1)/x_tilde.h(3),x_tilde.h(2)/x_tilde.h(3),'+k','MarkerSize',10,'LineWidth',1);

    xlim_aE=xlim;
    ylim_aE=ylim;

    subplot(3,2,3)
    hold on
    title('ALGs')
    for n = 1:N_samples
        plot(est_x_mean_a(1)+(est_x_samples_a(n,1)-est_x_mean_a(1)),...
            est_x_mean_a(2)+(est_x_samples_a(n,2)-est_x_mean_a(2)),'.k');
    end
    sugr_plot_Point_2D(est_p_a_emp,'ok','--k',1,3);

    plot(est_x_mean_a(1),est_x_mean_a(2),'ow','MarkerSize',12,'LineWidth',7);
    plot(est_x_mean_a(1),est_x_mean_a(2),'ok','MarkerSize',12,'LineWidth',3);

    plot(x_tilde.h(1)/x_tilde.h(3),x_tilde.h(2)/x_tilde.h(3),'+k','MarkerSize',10,'LineWidth',1);

    xlim_a = xlim;
    ylim_a = ylim;

    % geometric solution
    subplot(3,2,4)
    hold on
    title('SSDe')
    for n = 1:N_samples
        plot(est_x_mean_g(1)+(est_x_samples_g(n,1)-est_x_mean_g(1)),...
            est_x_mean_g(2)+(est_x_samples_g(n,2)-est_x_mean_g(2)),'.k');
    end

    sugr_plot_Point_2D(est_p_g_emp,'ok','--k',1,3);

    plot(est_x_mean_g(1),est_x_mean_g(2),'ow','MarkerSize',12,'LineWidth',7);
    plot(est_x_mean_g(1),est_x_mean_g(2),'ok','MarkerSize',12,'LineWidth',3);

    plot(x_tilde.h(1)/x_tilde.h(3),x_tilde.h(2)/x_tilde.h(3),'+k','MarkerSize',10,'LineWidth',1);

    xlim_g = xlim;
    ylim_g = ylim;

    % ml-estimation Hessian

    subplot(3,2,5)
    hold on
    title('MLEe')
    for n = 1:N_samples
        plot(est_x_mean_mH(1)+(est_x_samples_mH(n,1)-est_x_mean_mH(1)),...
            est_x_mean_mH(2)+(est_x_samples_mH(n,2)-est_x_mean_mH(2)),'.k');
    end
    sugr_plot_Point_2D(est_p_mH_emp,'ok','--k',1,3);
    plot(est_x_mean_mH(1),est_x_mean_mH(2),'ow','MarkerSize',12,'LineWidth',7);
    plot(est_x_mean_mH(1),est_x_mean_mH(2),'ok','MarkerSize',12,'LineWidth',3);

    plot(x_tilde.h(1)/x_tilde.h(3),x_tilde.h(2)/x_tilde.h(3),'+k','MarkerSize',10,'LineWidth',1);

    xlim_mH = xlim;
    ylim_mH = ylim;

    % ml-estimation        
    subplot(3,2,6)
    hold on
    title('MLEs')
    for n = 1:N_samples
        plot(est_x_mean_ml(1)+(est_x_samples_ml(n,1)-est_x_mean_ml(1)),...
            est_x_mean_ml(2)+(est_x_samples_ml(n,2)-est_x_mean_ml(2)),'.k');
    end
    sugr_plot_Point_2D(est_p_ml_emp,'ok','--k',1,3);

    plot(est_x_mean_ml(1),est_x_mean_ml(2),'ow','MarkerSize',12,'LineWidth',7);
    plot(est_x_mean_ml(1),est_x_mean_ml(2),'ok','MarkerSize',12,'LineWidth',3);

    plot(x_tilde.h(1)/x_tilde.h(3),x_tilde.h(2)/x_tilde.h(3),'+k','MarkerSize',10,'LineWidth',1);

    % determine common plot bounding box
    xlim_ml = xlim;
    ylim_ml = ylim;

    mxlimg = [min([xlim_a(1),xlim_g(1)]),...
        max([xlim_a(2),xlim_g(2)])];
    mylimg = [min([ylim_a(1),ylim_g(1)]),...
        max([ylim_a(2),ylim_g(2)])];

    subplot(3,2,2);xlim(mxlimg);ylim(mylimg);
    subplot(3,2,4);xlim(mxlimg);ylim(mylimg);

    mxliml = [min([xlim_mH(1),xlim_ml(1)]),...
        max([xlim_mH(2),xlim_ml(2)])];
    myliml = [min([ylim_mH(1),ylim_ml(1)]),...
        max([ylim_mH(2),ylim_ml(2)])];

    subplot(3,2,2);xlim(mxlimg);ylim(mylimg);
    subplot(3,2,4);xlim(mxlimg);ylim(mylimg);
    subplot(3,2,5);xlim(mxliml);ylim(myliml);
    subplot(3,2,6);xlim(mxliml);ylim(myliml);

    % same bb for all
    mx = [min([xlim_a(1),xlim_g(1),xlim_mH(1),xlim_ml(1)]),...
        max([xlim_a(2),xlim_g(2),xlim_mH(2),xlim_ml(2)])];
    my = [min([ylim_a(1),ylim_g(1),ylim_mH(1),ylim_ml(1)]),...
        max([ylim_a(2),ylim_g(2),ylim_mH(2),ylim_ml(2)]) ];

    subplot(3,2,2);xlim(mx);ylim(my);
    subplot(3,2,3);xlim(mx);ylim(my);
    subplot(3,2,4);xlim(mx);ylim(my);
    subplot(3,2,5);xlim(mx);ylim(my);
    subplot(3,2,6);xlim(mx);ylim(my);
    
    %% evaluate correctness
    S = 0.99;
    check_estimation_result(R,x_tilde_e,Cxx_a_e,ones(N_samples,1),est_x_samples_a,S,' ALGs');
    check_estimation_result(R,x_tilde_e,Cxx_mH_e,sigma_2_mH,est_x_samples_mH,S,' MLEe');
    set(gcf,'Position',[0.25*ss(1),70,0.35*ss(1),0.25*ss(2)])
    check_estimation_result(R,x_tilde_e,Cxx_ml_e,sigma_2_ml,est_x_samples_ml,S,' MLEs');
    set(gcf,'Position',[0.63*ss(1),70,0.35*ss(1),0.25*ss(2)])
    
    % Bias and  accuracy
    disp('----------------------------')
    disp('Precision, bias and accuracy')
    lb = max([trace(est_x_cov_aE),...
              trace(est_x_cov_a),...
              trace(est_x_cov_g),...
              trace(est_x_cov_ml),...
              trace(est_x_cov_mH)]);
    
    f = 10^-(ceil(log(lb)/log(10)/2));
    factor_for_bias=1/f;
    
    % provide bias, accuracy and test statistic
    disp('sigma = sqrt(tr(CovM)/2)')
    disp('bias  = true_x - est_x')
    disp('A     = |bias|^2 + sigma^2)')
    disp(['Values * ',num2str(f),':'])
    ALGe_bias_x________y_____sigma_________A = ...
        f*[-x_tilde_e'+est_x_mean_aE,sqrt(trace(est_x_cov_aE)/2),...
        sqrt(norm(x_tilde_e'-est_x_mean_aE)^2+trace(est_x_cov_aE)/2)]      %#ok<*NOPTS>
    
    ALGs_bias_x________y_____sigma_________A = ...
        f*[-x_tilde_e'+est_x_mean_a,sqrt(trace(est_x_cov_a)/2),...
        sqrt(norm(x_tilde_e'-est_x_mean_a)^2+trace(est_x_cov_a)/2)]
    
    SSDe_bias_x________y_____sigma_________A = ...
        f*[-x_tilde_e'+est_x_mean_g,sqrt(trace(est_x_cov_g)/2),...
        sqrt(norm(x_tilde_e'-est_x_mean_g)^2+trace(est_x_cov_g)/2)]
    
    MLEs_bias_x________y_____sigma_________A = ...
        f*[-x_tilde_e'+est_x_mean_mH,sqrt(trace(est_x_cov_mH)/2),...
        sqrt(norm(x_tilde_e'-est_x_mean_mH)^2+trace(est_x_cov_mH)/2)]
    
    MLEe_bias_x________y_____sigma_________A = ...
        f*[-x_tilde_e'+est_x_mean_ml,sqrt(trace(est_x_cov_ml)/2),...
        sqrt(norm(x_tilde_e'-est_x_mean_ml)^2+trace(est_x_cov_ml)/2)]
%%
end
% save('vp_0_200')