%% Fig. 4.11 test sensitivity factors in GMM with similarity
%
% see GMauss-Helmert model for similarity
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de

close all
clearvars
% clear all                                                                       %#ok<CLALL>

addpath(genpath('../../../General-Functions/'));
addpath(genpath('../Functions-GMM/'));

disp(' ')
disp(' ')
disp('===========r========================================================')
disp('---- Demo sensitivity factors in GMM with similarity Fig. 4.11 ----')
disp('-------------------------------------------------------------------')

%% Initiate parameters ####################################################

% init_rand = 0;
init_rand = 63;            % standard = 63
init_rand_seed(init_rand);

%type of generated data (1,2 or random) -----------------------------------

%data_type = -10;        % type of generated data
data_type =   2;         % standard = 2
if data_type == 0
    data_type = -10;
end

% statistical parameters --------------------------------------------------

% Significance level
S = 0.99;

% Tolerance
tol = chi2inv(S,2);

% selected parameters for sensititvity analysis ---------------------------
%
% 1,2 = translation (alternative choice)
% 3,4 = rotation, scale

rU = [3,4];

%% ########################################################################
disp('--------------------------------------------')
disp(' planar similarity with Gauss-Markov Model  ')
disp('--------------------------------------------')


%% generate data ##########################################################
% xis = a xi - b yi + c
% yis = b xi + a yi + d
switch data_type
    case 1
        tm  = [  1, 1; ...
            2, 2; ...
            -1,-2; ...
            -2,-1 ...
            ];
    case 2
        tm  = [  -7, 7 ;...
            1, 1; ...
            2, 2; ...
            -1,-2; ...
            -2,-1; ...
            ];
    otherwise
        tm=rand(-data_type,2)*5;
        
end
bb   = [min(tm(:,1)),max(tm(:,1)),min(tm(:,2)),max(tm(:,2))];
maxd = max(max(tm(:,1))-min(tm(:,1)),max(tm(:,2))-min(tm(:,2)));

% true parameters
xt  = [2,0.5,3,-2]';
disp('Parameters')
disp('   [ x(1)     -x(2)      x(3) ]');
disp('   [ x(2)      x(1)      x(4) ]');
True_transformation = [xt(1) -xt(2) xt(3);xt(2) xt(1) xt(4)] %#ok<*NOPTS>

%tm=rand(8,2)*2;
N   = length(tm(:));
d_I = 2;
I   = N/d_I;
Number_of_points = I
U   = 4;
sigma  = 0.01;
Standard_deviation_observations_m = sigma


%% setup elements of estimation
% observed coordinates
lm = zeros(I,d_I);
% design matrix A =[A_i']
Am = spalloc(N,U,N*U);
Cov_ll_m = zeros(I,I-1);
true_error = zeros(I,d_I);
for i =1:I
    true_error(i,:) = randn(2,1)*sigma;
    lm(i,:)         = [xt(1) -xt(2) xt(3);xt(2) xt(1) xt(4)]*[tm(i,:)';1] ...
        + true_error(i,:)';
    Am(2*i-1:2*i,:) = [tm(i,1) -tm(i,2) 1 0;...
        tm(i,2)  tm(i,1) 0 1];                               %#ok<*SPRIX>
    Cov_ll_m(i,:)   = [1 0 0 1]*sigma^2;
end
true_errors=true_error
true_and_observed_coordinates=[tm , lm]
av  = zeros(N,1);

%% factors for plotting
factor_v  = 0.5*maxd/sqrt(N)/sigma;
factor_r  = 0.5*maxd/sqrt(N);
factor_X  = maxd/sqrt(N)/6;
factor_mu = maxd/sqrt(N)/4;
factor_nv = factor_v/10;

%% estimate parameters

[est_x,Cov_xx,sigma_0q,R,vm,Xv,Rim,nabla_lv,muv,muv1,Um1q,Um2]=...
    GaussMarkovModelLinear_groups(lm,Cov_ll_m,Am,av,rU);



%% plot results
screensize = plot_init;
figure('Color','w','Position',[100 100 screensize+[-300 -300]]);

subplot(2,3,1);
hold on;
plot(tm(:,1),tm(:,2),'.r','Markersize',12);       % given data
for i=1:I
    plot([tm(i,1),tm(i,1)+factor_v*vm(i,1)],[tm(i,2),tm(i,2)+factor_v*vm(i,2)],'-k');
end
title({'residuals [m]', strcat('max=',num2str(max(norm(vm(:,1:2)))))});
xlim([min(tm(:,1))-1,max(tm(:,1))+1]);
ylim([min(tm(:,2))-1,max(tm(:,2))+1]);
axis equal;


subplot(2,3,2);
hold on;
for i=1:I
    plot_circle(tm(i,1),tm(i,2),factor_r*Rim(i,1),'-r');
end
title({'redundancy number', strcat('min=',num2str(min(Rim(:,1))))});
xlim([min(tm(:,1))-factor_r*1,max(tm(:,1))+factor_r*1]);
ylim([min(tm(:,2))-factor_r*1,max(tm(:,2))+factor_r*1]);
axis equal;


subplot(2,3,4);
hold on;
for i=1:I
    plot_circle(tm(i,1),tm(i,2),factor_X*sqrt(Xv(i)/d_I),'-k','LineWidth',2);
    plot_circle(tm(i,1),tm(i,2),factor_X*sqrt(tol/d_I),'--r');
end
title({'test statistics', strcat('max=',num2str(max(sqrt(Xv/d_I)))),'- - - critical values'});
xlim([min(tm(:,1))-factor_X*5,max(tm(:,1))+factor_X*5]);
ylim([min(tm(:,2))-factor_X*5,max(tm(:,2))+factor_X*5]);
axis equal;

subplot(2,3,3);
hold on;
for i=1:I
    plot_circle(tm(i,1),tm(i,2),factor_mu*muv(i),'-k');
end
title({'sensitivity factors', strcat('max=',num2str(max(muv)))});
xlim([min(tm(:,1))-factor_mu*5,max(tm(:,1))+factor_mu*5]);
ylim([min(tm(:,2))-factor_mu*5,max(tm(:,2))+factor_mu*5]);
axis equal;


subplot(2,3,5);
hold on;
for i=1:I
    plot_circle(tm(i,1),tm(i,2),factor_nv*nabla_lv(i),'-r');
end
title({'min. detectable outliers [m]';strcat('max=',num2str(max(nabla_lv)))});
xlim([min(tm(:,1))-1,max(tm(:,1))+1]);
ylim([min(tm(:,2))-1,max(tm(:,2))+1]);
axis equal;

if ~isempty(rU)
    subplot(2,3,6);
    hold on;
    for i=1:I
        plot_circle(tm(i,1),tm(i,2),factor_mu*muv1(i),'-k');
    end
    title({'sens. factors partial';strcat('max=',num2str(max(muv1)),', param=',num2str(rU))});
    xlim([min(tm(:,1))-factor_mu*5,max(tm(:,1))+factor_mu*5]);
    ylim([min(tm(:,2))-factor_mu*5,max(tm(:,2))+factor_mu*5]);
    axis equal;
end


%%

disp('.............................')
disp('        diagnostics          ')
disp('.............................')
Estimated_transformation = ...
    [est_x(1) -est_x(2) est_x(3);est_x(2) est_x(1) est_x(4)]
Theoretical_precision_of_paramters   = full(sqrt(diag(Cov_xx)))'
Redundancy              =R
Estimated_sigma_0       = sqrt(sigma_0q)

a________I________vx________vy_________z=...
    [(1:I)',vm,sqrt(Xv/d_I)]

[m,i]=max(sqrt(diag(true_error*true_error')));
disp(horzcat('Maximal true error                 = ',num2str(m,'% 12.5f'),...
    ' at: ', num2str(i,'% 3.0f')));
[m,i]=max(sqrt(diag(vm*vm')));
disp(horzcat('Maximal residual                   = ',num2str(m,'% 12.5f'),...
    ' at: ', num2str(i,'% 3.0f')));
[m,i]=min(Rim(:,1));
disp(horzcat('Minimal redundancy number          = ',num2str(m,'% 12.5f'),...
    ' at: ', num2str(i,'% 3.0f')));

[m,i]=max(sqrt(Xv/d_I));
if max(Xv) > tol
    disp(horzcat('Maximal test statistic             = ',num2str(m,'% 12.5f'),...
        ' at: ', num2str(i,'% 3.0f')), ' ***');
else
    disp(horzcat('Maximal test statistic             = ',num2str(m,'% 12.5f'),...
        ' at: ', num2str(i,'% 3.0f')));
end
[m,i]=max(nabla_lv);
disp(horzcat('Max. of minimal detectable outlier = ',num2str(m,'% 12.5f'),...
    ' at: ', num2str(i,'% 3.0f')));

[m,i]=max(muv);
disp(horzcat('Max. sensitivity factor            = ',num2str(m,'% 12.5f'),...
    ' at: ', num2str(i,'% 3.0f')));
if ~isempty(rU)
    [m,i]=max(muv1);
    disp(horzcat('Max. sensitivity factor [',num2str(rU),']     = ',num2str(m,'% 12.5f'),...
        ' at: ', num2str(i,'% 3.0f')));
end



