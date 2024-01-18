%% test whether selection is necessary when using L through X constraint
% Result: they differ, but not statistically significant
%
% Case o: no selection
% Case m: with selection
%
% Method
% 1. Generate true values for X and three L_i = X \wedge Y_i
% 2. Derive theoretical covariance matrix for both cases
% 3. Generate N_iter samples with random perturbations
% 4. Determine solution for both cases
% 5. determine empirical covariance matrix for both cases
% 6. Statistically test empirical vs. theoretical covariance matrix
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 

% clear all                                                                  %#ok<CLALL>
clearvars
close all
clc

addpath(genpath('../General-Functions'))

S_level = 0.99; % significance level
sigma = 0.001;  % noise standard deviation for perturbation of X
N_iter = 500;   % number of samples, can be enlarged for stronger testing

disp('----- Test: selection of constraints for estimating X from (Li) -----')

X_true = randn(4,1);             
X_true = X_true/norm(X_true);         % spherically normalized X

Y1 = randn(4,1);                      
Y2 = randn(4,1);
Y3 = randn(4,1);                      % three other points Y_i

L1_true = calc_Pi(X_true)*Y1;         
L2_true = calc_Pi(X_true)*Y2;
L3_true = calc_Pi(X_true)*Y3;
L1_true = L1_true/norm(L1_true);
L2_true = L2_true/norm(L2_true);
L3_true = L3_true/norm(L3_true);      % spherically normalized L_i

Cov1 = sigma^2*(eye(6)-...
    [L1_true,calc_Dual*L1_true]*[L1_true,calc_Dual*L1_true]');
Cov2 = sigma^2*(eye(6)-...
    [L2_true,calc_Dual*L2_true]*[L2_true,calc_Dual*L2_true]');
Cov3 = sigma^2*(eye(6)-...
    [L3_true,calc_Dual*L3_true]*[L3_true,calc_Dual*L3_true]');
Cov = [Cov1 zeros(6,12);...
     zeros(6,6) Cov2 zeros(6,6);...
     zeros(6,12) Cov3];             % covariance matrix of L_i
%% Case without selection; determine theoretical CovM
Ao = [calc_Gammadual(L1_true);...
    calc_Gammadual(L2_true);...
    calc_Gammadual(L3_true)];       % Coefficient matrix 

[Uo,Do,Vo] = svd(Ao);                 % SVD
 
Apo = Vo(:,1:3)*inv(Do(1:3,1:3))*Uo(:,1:3)';                                 %#ok<*MINV>
                                    % Pseudo inverse
                                    
B = [calc_Pidual(X_true)' zeros(4,6) zeros(4,6);...
     zeros(4,6) calc_Pidual(X_true)' zeros(4,6);...
     zeros(4,6) zeros(4,6) calc_Pidual(X_true)'];
                                    % Jacobian
CovMo = Apo*(B*Cov*B')*Apo';          % Theoretical CovM of estX

Jr = null(X_true');
CovM_rr_o = Jr'*CovMo*Jr;            % Theoretical CovM of estX_r
disp('Theoretical covariance matrix of estX_r, without selection: ');disp(CovM_rr_o)

%% Case with selection; determine theoretical CovM
[Gr1,M1] = calc_Gammadual_reduced(L1_true);
[Gr2,M2] = calc_Gammadual_reduced(L2_true);
[Gr3,M3] = calc_Gammadual_reduced(L3_true);
Am = [Gr1;Gr2;Gr3];                   % Coefficient matrix 

[Um,Dm,Vm] = svd(Am);                 % SVD

Apm = Vm(:,1:3)*inv(Dm(1:3,1:3))*Um(:,1:3)';
                                    % Pseudo inverse
                                    
B = [M1*calc_Pidual(X_true)' zeros(2,6) zeros(2,6);...
     zeros(2,6) M2*calc_Pidual(X_true)' zeros(2,6);...
     zeros(2,6) zeros(2,6) M3*calc_Pidual(X_true)'];
                                    % Jacobian
 
CovMm = Apm*(B*Cov*B')*Apm';          % Theoretical CovM of estX
CovM_rr_m = Jr'*CovMm*Jr;            % Theoretical CovM of estX_r
disp('Theoretical covariance matrix of estX_r, with selection: ');disp(CovM_rr_m)



%% Simulate N_iter cases
mean_o = zeros(3,1);
mean_m = zeros(3,1);
CovM_emp_o = zeros(3);
CovM_emp_m = zeros(3);                % initiate empirical CovM

for iter = 1:N_iter
    L1 = rand_gauss(L1_true,Cov1,1);
    L2 = rand_gauss(L2_true,Cov2,1);
    L3 = rand_gauss(L3_true,Cov3,1);  % add noise
    L1 = sugr_constrain_3D_Line(L1);
    L2 = sugr_constrain_3D_Line(L2);
    L3 = sugr_constrain_3D_Line(L3);% enforce Plücker constraint
    L1 = L1/norm(L1);
    L2 = L2/norm(L2);
    L3 = L3/norm(L3);                 % normalize spherically         

    % without selection
    A = [calc_Gammadual(L1);calc_Gammadual(L2);calc_Gammadual(L3)];
                                    % Coefficient matrix 
                                    
    [~,~,V] = svd(A);                 % SVD
    
    X_est_o = V(:,4);                % algebraically best point
    
    Xr_o = Jr' * X_est_o;             % reduced coordinates
    
    mean_o = mean_o+Xr_o;
    CovM_emp_o = CovM_emp_o+Xr_o*Xr_o';% update empirical CovM
    
    % with selection, use same simulated data
    Gr1 = calc_Gammadual_reduced(L1);
    Gr2 = calc_Gammadual_reduced(L2);
    Gr3 = calc_Gammadual_reduced(L3);
    A = [Gr1;Gr2;Gr3];                
                                    % Coefficient matrix 
    [U,D,V] = svd(A);                 % SVD
    
    X_est_m = V(:,4);                % algebraically best X
    Xr_m = Jr' * X_est_m;             % reduced coordinates
    mean_m = mean_m+Xr_m;
    CovM_emp_m = CovM_emp_m+Xr_m*Xr_m';% update empirical CovM
end
mean_o = mean_o/N_iter;
mean_m = mean_m/N_iter;
CovM_emp_m = CovM_emp_m/N_iter;
CovM_emp_o = CovM_emp_o/N_iter;        % normalize CovM

disp('Empirical covariance matrix of estX_r, with selection: ');disp(CovM_emp_m)
disp('Empirical covariance matrix of estX_r, without selection: ');disp(CovM_emp_o)

disp(strcat('Sample size: ',num2str(N_iter)));
disp('Test of bias')
Xo = mean_o'*inv(CovM_rr_o)*mean_o*N_iter;
Xm = mean_m'*inv(CovM_rr_m)*mean_m*N_iter;
Tm = chi2inv(S_level,3);
To = Tm;
if Xo < To
    disp(strcat('Estimation without selection : ',num2str(Xo),' < [',num2str(To),']'));
else  
    disp(strcat('Estimation without selection : ',num2str(Xo),' > [',num2str(To),'] ***'));
end
if Xm < Tm
    disp(strcat('estimation with selection    : ',num2str(Xm),' < [',num2str(Tm),']'));
else  
    disp(strcat('estimation with selection    : ',num2str(Xm),' > [',num2str(Tm),'] ***'));
end
disp(' ')
disp('Test of empirical and theoretical CovM')
[lambdam,Tm] = check_CovM(CovM_emp_m,CovM_rr_m,N_iter,S_level);
[lambdao,To] = check_CovM(CovM_emp_o,CovM_rr_o,N_iter,S_level);
                                    % test statistic and F-fractile
lo_To_lm_Tm = [lambdao,lambdam,Tm];
if lambdao > To(1) && lambdao < To(2)
    disp(strcat('estimation without selection : ',num2str(lambdao),'     in [',num2str(To),']'));
else  
    disp(strcat('estimation without selection : ',num2str(lambdao),' not in [',num2str(To),'] ***'));
end
if lambdam > To(1) && lambdao < To(2)
    disp(strcat('estimation with selection    : ',num2str(lambdam),'     in [',num2str(Tm),']'));
else  
    disp(strcat('estimation with selection    : ',num2str(lambdam),' not in [',num2str(Tm),'] ***'));
end




