%% GHM with similarity
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% see GMM for similarity
% See also demo_GHM_2D_lines

close all
clearvars
% clear all                                                                  %#ok<CLALL>

addpath(genpath('../../../General-Functions/'));
addpath(genpath('../Functions-GHM/'));
 
simulate_GMM_similarity = true;  % choose transformed coordinates CovM=0
%simulate_GMM_similarity = false; % choose transformed coordinates CovM~=0

 % initialization for random numbers
init_rand = 63;            
init_rand_seed(init_rand);

%data_type = -10;         % type of generated data
data_type =   2;         % type of generated data
if data_type == 0
    data_type = -10;
end

% statistical parameters
S = 0.99;                % Significance level
tol = chi2inv(S,2);

% controlling convergence
Tx      = 0.00001; % threshold for convergence
maxiter = 10;   % maximum number of iterations

% selected parameters for sensititvity analysis
%
% 1,2 = translation (alternative choice)
% 3,4 = rotation, scale
rU = [3,4];              


disp('============================================')
disp(' Planar similarity with Gauss-Helmert model ')
disp('--------------------------------------------')

if simulate_GMM_similarity
    disp('Simulate similarity with Gauss-Markov model');
    disp('*******************************************')
end

%% model ------------------------------------------------------------------
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
% for plotting: bounding box
bb   = [min(tm(:,1)),max(tm(:,1)),min(tm(:,2)),max(tm(:,2))];
maxd = max(max(tm(:,1))-min(tm(:,1)),max(tm(:,2))-min(tm(:,2)));

xt  = [2,0.5,3,-2]';
disp('Parameters')
disp('   [ x(1)     -x(2)      x(3) ]');
disp('   [ x(2)      x(1)      x(4) ]');
True_transformation = [xt(1) -xt(2) xt(3);xt(2) xt(1) xt(4)]                %#ok<NOPTS>

% number of observations
N   = 2*length(tm(:));
% number of observations per correspondence
d_I = 4;
% number of correspondences = numbe rof points
I   = N/d_I;
Number_of_points = I                                                        %#ok<NOPTS>
% number of unknown parameters
U   = 4;
% standard deviations for generation and estimation
sigma_x = 0.01;                                                             %#ok<NASGU>
sigma_y  = 0.01;
if simulate_GMM_similarity
    sigma_x  = 0.00;
end
Standard_deviation_observations_m = [sigma_x,sigma_y]                       %#ok<NOPTS>

%% generate observations --------------------------------------------------
lm         = zeros(I,d_I);           % Ix4 matrix
true_error = zeros(I,d_I);
Cov_ll_m   = zeros(I,d_I^2);         % I x 16 matrix
for i =1:I
    % random perturbations
    if simulate_GMM_similarity
         true_error(i,:)  = [[0 0],randn(1,2)*sigma_y]; 
    else
        true_error(i,:) = [randn(1,2)*sigma_x,randn(1,2)*sigma_y];          %#ok<UNRCH>
    end
    
    lm(i,:)= [tm(i,:),...
             ([xt(1) -xt(2) xt(3);xt(2) xt(1) xt(4)]*[tm(i,:)';1])'] ...
           + true_error(i,:);
    Am(2*i-1:2*i,:) = [tm(i,1) -tm(i,2) 1 0;...
                       tm(i,2)  tm(i,1) 0 1];
    Cov_ll_m(i,:)= [sigma_x^2*[1 0 0 0 0 1 0 0], ...
                    sigma_y^2*[0 0 1 0 0 0 0 1]];
end
true_errors = true_error                                                    %#ok<NOPTS>
observed_coordinates = lm                                                   %#ok<NOPTS>

%% factors for plotting ---------------------------------------------------
factor_v  = 0.5*maxd/sqrt(N/2)/max([sigma_x,sigma_y]);
factor_r  = 0.5*maxd/sqrt(N/2);
factor_X  = maxd/sqrt(N/2)/6;
factor_mu = maxd/sqrt(N/2)/4;
factor_nv = factor_v/20;

%% estimate parameters ----------------------------------------------------

% xa = [1.8,0.6,2.5,-3]';    % approximate values, true [2,0.5,3,-2]';
xa = zeros(4,1);
sa = max([sigma_x,sigma_y]) * [ [1,1]/maxd, [1,1]/sqrt(I) ];

[est_x,Cov_xx,sigma_0q,R,vm,Xv,Rim,nabla_lv,muv,muv1,uv1q,uv2] = ...
            GaussHelmertModelGroups(lm, Cov_ll_m, @cg_similarity_2D, xa, ...
            @ux_similarity_2D, sa, Tx, maxiter, rU);
          
%% plot results ----------------------------------------------------------- 

% for normalizing the test statistic of the observations
% assuming it ist the same for all observations
rk_Cov_ll = rank(reshape(Cov_ll_m(1,:),4,4));

ScrS = plot_init;  
figure('Color','w','Position',[100 100  ScrS+[ -600 -300]])
subplot(2,3,1)
hold on
plot(tm(:,1),tm(:,2),'.r','Markersize',12)       % given data
for i = 1:I
    plot([tm(i,1),tm(i,1)+factor_v*(vm(i,3))],[tm(i,2),tm(i,2)+factor_v*(vm(i,4))],'-k');
end
title({'residuals [m]', strcat('max=',num2str(max(norm(vm(:,1:4)))))})
xlim([min(tm(:,1))-1,max(tm(:,1))+1]);
ylim([min(tm(:,2))-1,max(tm(:,2))+1]);
axis equal


subplot(2,3,2)
hold on
for i = 1:I
    plot_circle(tm(i,1),tm(i,2),factor_r*Rim(i,1),'-r');
end
title({'redundancy number', strcat('min=',num2str(min(Rim(:,1))))})
xlim([min(tm(:,1))-factor_r*1,max(tm(:,1))+factor_r*1]);
ylim([min(tm(:,2))-factor_r*1,max(tm(:,2))+factor_r*1]);
axis equal

subplot(2,3,4)
hold on
for i=1:I
    plot_circle(tm(i,1),tm(i,2),factor_X*sqrt(Xv(i)/rk_Cov_ll),'-k','LineWidth',2);
    plot_circle(tm(i,1),tm(i,2),factor_X*sqrt(tol/rk_Cov_ll),'--r');
end
title({'test statistics', strcat('max=',num2str(max(sqrt(Xv/rk_Cov_ll)))),'- - - critical values'});
xlim([min(tm(:,1))-factor_X*5,max(tm(:,1))+factor_X*5]);
ylim([min(tm(:,2))-factor_X*5,max(tm(:,2))+factor_X*5]);
axis equal

subplot(2,3,3)
hold on
for i = 1:I
    plot_circle(tm(i,1),tm(i,2),factor_mu*muv(i),'-k');
end
title({'sensitivity factors', strcat('max=',num2str(max(muv)))})
xlim([min(tm(:,1))-factor_mu*5,max(tm(:,1))+factor_mu*5]);
ylim([min(tm(:,2))-factor_mu*5,max(tm(:,2))+factor_mu*5]);
axis equal

subplot(2,3,5)
hold on
for i = 1:I
    plot_circle(tm(i,1),tm(i,2),factor_nv*nabla_lv(i),'-r');
end
title({'min. detectable outliers [m]';strcat('max=',num2str(max(nabla_lv)))})
xlim([min(tm(:,1))-1,max(tm(:,1))+1]);
ylim([min(tm(:,2))-1,max(tm(:,2))+1]);
axis equal

if ~isempty(rU)
    subplot(2,3,6)
    hold on
    for i = 1:I
        plot_circle(tm(i,1),tm(i,2),factor_mu*muv1(i),'-k');
    end
    title({'sens. factors partial';strcat('max=',num2str(max(muv1)),', param=',num2str(rU))})
    xlim([min(tm(:,1))-factor_mu*5,max(tm(:,1))+factor_mu*5]);
    ylim([min(tm(:,2))-factor_mu*5,max(tm(:,2))+factor_mu*5]);
    axis equal
end

%% diagnostics ------------------------------------------------------------
disp('diagnostics                  ')
disp('.............................')
Estimated_transformation = ...
    [est_x(1) -est_x(2) est_x(3);est_x(2) est_x(1) est_x(4)]                %#ok<NOPTS>
Theoretical_precision   = sqrt(diag(full(Cov_xx)))'                         %#ok<NOPTS>
Redundancy              = R                                                 %#ok<NOPTS>
Estimated_sigma_0       = sqrt(sigma_0q)                                    %#ok<NOPTS>
a________I_______vxl_______vyl_______vxr_______vyr_________X = ...
[(1:I)',vm,sqrt(Xv/rk_Cov_ll)]                                              %#ok<NOPTS>

[m,i] = max(sqrt(diag(true_error*true_error')));
    disp(horzcat('Maximal true error                 = ',num2str(m,'% 12.5f'),...
                   ' at: ', num2str(i,'% 3.0f')));
[m,i] = max(sqrt(diag(vm*vm')));
    disp(horzcat('Maximal residual                   = ',num2str(m,'% 12.5f'),...
                   ' at: ', num2str(i,'% 3.0f')));
[m,i] = min(Rim(:,1));
    disp(['Minimal redundancy number          = ',num2str(m,'% 12.5f'),...
                   ' at: ', num2str(i,'% 3.0f')]);

[m,i] = max(sqrt(Xv/rk_Cov_ll));
if max(Xv) > tol 
    disp(['Maximal test statistic             = ',num2str(m,'% 12.5f'),...
                   ' at: ', num2str(i,'% 3.0f'), ' ***']);
else
    disp(['Maximal test statistic             = ',num2str(m,'% 12.5f'),...
                   ' at: ', num2str(i,'% 3.0f')]);
end
[m,i] = max(nabla_lv);
    disp(['Max. of minimal detectable outlier = ',num2str(m,'% 12.5f'),...
                   ' at: ', num2str(i,'% 3.0f')]);
               
[m,i] = max(muv);  
    disp(['Max. sensitivity factor            = ',num2str(m,'% 12.5f'),...
                   ' at: ', num2str(i,'% 3.0f')]);
if ~isempty(rU) 
    [m,i] = max(muv1); 
    disp(['Max. sensitivity factor [',num2str(rU),']     = ',num2str(m,'% 12.5f'),...
                   ' at: ', num2str(i,'% 3.0f')]);
end



