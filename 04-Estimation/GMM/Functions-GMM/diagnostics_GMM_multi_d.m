%% Diagnostics Gauss-Markov model multidimensional
%
% given all parameters of a linearized estimation problem
% determine the essential diagnostic values
% assuming the number N of observations is small enough
%    (requires several NxN matrices)
%
% r_U    = index set for parameters to be estimated (others are nuisance)
% I      = number of observational groups
% d_I    = dimension of observatioanl groups (must be the same for all)
% Am     = Jacobian
% Cov_ll = CovM of observations
% W_ll   = WeightM of observations
% Cov_xx = CovM of all estimated parameters
% W_xx   = WeightM of all estimated parameters
% vv     = residuals
%
%
% Xv     = standardized test statistics (chi^2)
% Rim    = redundancy matrices as row vectors (I * d_I^2)
% nabla_lv = lower bound for detectable outliers (maximum eigenvalue)
% muv      = sensitivity factor for all parameters
% muv1     = sensitivity factor for parameters specified by r_U
% Um1q     = diagonal blocks of \bar U_1 (numerator of muv1^2)
% Um2      = diagonal blocks of U2
%
%
% Wolfgang Förstner 10/2016
% wfoerstn@uni-bonn.de


function [vi,Xv,Rim,nabla_lv,muv,muv1,U1mq,U2m] = diagnostics_GMM_multi_d...
    (r_U,I,d_I,Am,Cov_ll,W_ll,Cov_xx,W_xx,vv)

[~,U] = size(Am);

delta_0 = 4.13;  % all sizes referring to d=1.

if ~isempty(r_U)
    % ranges
    r_C = r_U;
    r_D = setdiff(1:U,r_U);
    % partitioned design matrix A=[C,D]
    Cm = Am(:,r_C);
    Dm = Am(:,r_D);
    % reduced design matrix
    Cmr    = Cm - Dm * inv(W_xx(r_D,r_D)) * W_xx(r_D,r_C);  % C_reduced
    Cov_11 = Cov_xx(r_C,r_C);
end
%% observational groups
Rim  = zeros(I,d_I^2);
U1mq = zeros(I,d_I^2);
U2m  = zeros(I,d_I^2);
vi   = zeros(I,d_I);
Xv   = zeros(I,1);
nabla_lv = zeros(I,1);
muv      = zeros(I,1);
muv1     = zeros(I,1);
for i=1:I
    % detectability, test statistics, sensitivity wrt all parameters
    i_range     = d_I*i-(d_I-1):d_I*i;
    Ai          = Am(i_range,:)';
    Cov_ll_ii   = Cov_ll(i_range,i_range);  % CovM of group
    W_ll_ii     = W_ll(i_range,i_range);    % Weight matrix of group
    Cov_vv_ii   = Cov_ll_ii-Ai'*Cov_xx*Ai;  % covariance matrix of resiudal
    W_vv_ii     = inv(Cov_vv_ii+eps);       % Weihgt matrix of residuals
    Rii         = Cov_vv_ii*W_ll_ii;        % redundancy matrix
    Rim(i,:)    = Rii(:)';                  % vectorized Rii for storing
    Rii_inv     = full(inv(Rii));           % inverse of redundancy matrix
    vi(i,:)     = vv(i_range)';             % residual of group
    Xv(i)       = vv(i_range)'*W_vv_ii* vv(i_range);      %#ok<*MINV>
    % Xhi^2-test statistic
    nabla_lv(i) = delta_0 * ...
        sqrt(real(max(eig(Rii_inv*Cov_ll_ii))));
    % minimum detectable outlier
    muv(i)      = sqrt(real(max(eig(Rii_inv-eye(d_I)))));
    % sensitivity factors
    
    % sensitivity wrt selected parameters range rU
    if ~isempty(r_U)
        U1q       = Cmr(i_range,:) * Cov_11  * Cmr(i_range,:)' * W_ll_ii;
        U1mq(i,:) = U1q(:)';
        U2        = Dm(i_range,:)  / W_xx(r_D,r_D)  * ...
            Dm(i_range,:)' * W_ll_ii;
        U2m(i,:)  = U2(:)';
        %check     = U1q + U2 + Rii-eye(d_I)         % I=U_C_bar+U_D+R
        %i,check-eye(d_I)
        muv1(i)   = sqrt(real(max(eig(U1q * Rii_inv))));
        % sensitivity factors
    else
        U1mq=0;
        U2m =0;
        muv1=0;
    end
end
