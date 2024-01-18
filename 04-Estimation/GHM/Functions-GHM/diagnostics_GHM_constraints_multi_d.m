%% Diagnostics Gauss-Helmert model multidimensional
%
% referring to constraints
% assuming each group of observations belongs to one group of constraints
%      number = I, with d_I observations and d_G constraints
%
%  [vgi,Xv,Rgm,nabla_cg,muv,muv1,U1mqc,U2mc]=
%            diagnostics_GHM_constraints_multi_d...
%   (r_U,I,d_I,d_G,Am,Bm,Cov_xx,W_xx,W_gg,vv)
% 
% r_U    = index set for parameters to be estimated (others are nuisance)
% I      = number of observational groups
% d_I    = dimension of observatioanl groups (must be the same for all)
% Am     = G x U Jacobian
% Bm     = G x N Jacobain
% Cov_xx = U x U CovM of all estimated parameters
% W_xx   = U x U WeightM of all estimated parameters
% Wgg    = 
% vv     = N x 1 residuals
%
% vgi      = I x dI    matrix of residuals
% Xv       = I x 1     vector standardized test statistics (chi^2)
% Rim      = I * d_I^2 redundancy matrices as row vectors 
% nabla_lv = I x 1     lower bound for detectable outliers (maximum eigenvalue)
% muv      = I x 1     sensitivity factor for all parameters
% muv1     = I x 1     sensitivity factor for parameters specified by r_U
% Um1q     = I * d_I^2 diagonal blocks of \bar U_1 (numerator of muv1^2)
% Um2      = I * d_I^2 diagonal blocks of U2
%
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
% 
% See also GaussHelmertModel

function [vgi,Xv,Rgm,nabla_cg,muv,muv1,U1mqc,U2mc]=diagnostics_GHM_constraints_multi_d...
    (r_U,I,d_I,d_G,Am,Bm,Cov_xx,W_xx,W_gg,vv)

[~,U] = size(Am);

delta_0=4.13;  % all sizes referring to d_J.

if ~isempty(r_U)
    % ranges
    r_C= r_U;
    r_D= setdiff(1:U,r_U);
    % partitioned design matrix A=[C,D]
    Cm = Am(:,r_C);
    Dm = Am(:,r_D);
    % reduced design matrix
    Cmr    = Cm - Dm / W_xx(r_D,r_D) * W_xx(r_D,r_C);  % C_reduced
    Cov_11 = Cov_xx(r_C,r_C);
end

%% observational groups
Rgm   = zeros(I,d_G^2);
U1mqc = zeros(I,d_G^2);
U2mc  = zeros(I,d_G^2);
Xv    = zeros(I,1);
vgi   = zeros(I,d_G);
nabla_cg = zeros(I,1);
muv   = zeros(I,1);
muv1  = zeros(I,1);

for i=1:I
    % detectability, test statistics, sensitivity wrt all parameters
    i_range     = d_I*i-(d_I-1):d_I*i;
    g_range     = d_G*i-(d_G-1):d_G*i;
    Amt_i       = Am(g_range,:);
    Bmt_i       = Bm(i_range,g_range)';
    %W_gg_ii     = inv(Bmt_i * Cov_ll_ii * Bmt_i'); % Weight matrix of group
    W_gg_ii     = W_gg(g_range,g_range); % Weight matrix of group
    Cov_gg_ii   = inv(W_gg_ii);
    Cov_vgvg_ii = Cov_gg_ii - Amt_i * Cov_xx * Amt_i';
    W_vgvg_ii   = inv(Cov_vgvg_ii+eps);
    % redundancy matrix of constraints
    R_gg_ii     = full(Cov_vgvg_ii * W_gg_ii);
    R_gg_ii_inv = inv(R_gg_ii);
    Rgm(i,:)    = R_gg_ii(:)';                  % vectorized Rgg
    vgi(i,:)    = Bmt_i*vv(i_range);       % residual of group
    Xv(i)       = vgi(i,:) * W_vgvg_ii * vgi(i,:)';                         %#ok<MINV>
    % Xhi^2-test statistic
    nabla_cg(i) = delta_0 * ...
        sqrt(real(max(eig( R_gg_ii \ inv(W_gg_ii) ))));
    % minimum detectable outlier
    muv(i)      = sqrt(real(max(eig(R_gg_ii_inv - eye(d_G) ))));
    % sensitivity factors
    
    % sensitivity wrt selected parameters range rU
    
    if ~isempty(r_U)
        % U's for constraints
        U1qc       = Cmr(g_range,:) * Cov_11  * Cmr(g_range,:)' * W_gg_ii;
        U1mqc(i,:) = U1qc(:)';
        U2c        = Dm(g_range,:)  / W_xx(r_D,r_D)  * ...
            Dm(g_range,:)' * W_gg_ii;
        U2mc(i,:)  = U2c(:)';
        % check      = U1qc + U2c + R_gg_ii - eye(d_G)
        % I=U_C_bar+U_D+R, check-eye(d_I)
        
        % sensitivity factors
        muv1(i)   = sqrt(real(max(eig(U1qc / R_gg_ii ))));
    end
end

