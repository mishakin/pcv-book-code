%% Diagnostics Gauss-Helmert model multidimensional
%
% using sparse inverse
%
% params see GaussHelmertModel
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
% 
% See also GaussHelmertModel

function [vi,Xv,Rim,nabla_lv,muv,muv1,U1mq,U2m] = diagnostics_GHM_multi_d...
    (r_U,I,d_I,d_G,Am,Bm,Cov_ll,W_ll,Cov_xx,W_xx,vv)

[~,U] = size(Am);

delta_0=4.13;  % all sizes referring to d_I.

if ~isempty(r_U)
    % ranges
    r_C = r_U;
    r_D = setdiff(1:U,r_U);
    % partitioned design matrix A=[C,D]
    Cm  = Am(:,r_C);
    Dm  = Am(:,r_D);
    % reduced design matrix
    Cmr    = Cm - Dm / W_xx(r_D,r_D) * W_xx(r_D,r_C);  % C_reduced
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
    g_range     = d_G*i-(d_G-1):d_G*i;
    Am_i        = Am(g_range,:)';
    Bm_i        = Bm(i_range,g_range);
    Cov_ll_ii   = reshape(Cov_ll(i,:),d_I,d_I);  % CovM of group
    W_ll_ii     = W_ll(i_range,i_range);         % Weight matrix of group
    W_gg_ii     = Bm_i' * Cov_ll_ii * Bm_i;      % Weight matrix of group
    % Cov_ll*Bm*W_gg*Bm'-Cov_ll*Bm*W_gg*Am*Cov_xx*Am'*W_gg*Bm'
    Cov_vv_ii   = Cov_ll_ii * Bm_i *...
        (eye(d_G) - W_gg_ii *Am_i'*Cov_xx*Am_i)*...
        W_gg_ii*Bm_i'*Cov_ll_ii;  % covariance matrix of resiudal
    W_vv_ii     = inv(Cov_vv_ii+eps);       % Weight matrix of residuals
    Rii         = Cov_vv_ii*W_ll_ii;        % redundancy matrix
    Rim(i,:)    = Rii(:)';                  % vectorized Rii
    Rii_inv     = full(inv(Rii));           % inverse of redundancy matrix
    vi(i,:)     = vv(i_range)';             % residual of group
    Xv(i)       = vv(i_range)'*W_vv_ii* vv(i_range);      %#ok<MINV>
    % Xhi^2-test statistic
    nabla_lv(i) = delta_0 * ...
        sqrt(real(max(eig(Rii_inv*Cov_ll_ii))));
    % minimum detectable outlier
    muv(i)      = sqrt(real(max(eig(Rii_inv-eye(d_I)))));
    % sensitivity factors
    
    % sensitivity wrt selected parameters range rU
    if ~isempty(r_U)
        % U's for constraints
        U1qc       = Cmr(g_range,:) * Cov_11  * Cmr(g_range,:)' * W_gg_ii;
        U2c        = Dm(g_range,:)  / (W_xx(r_D,r_D))  * ...
            Dm(g_range,:)' * W_gg_ii;
        % check      = U1qc + U2c + (eye(G)-Am_i'*Cov_xx*Am_i*inv(Cov_gg_ii));
        % I=U_C_bar+U_D+R, check-eye(d_I)
        
        % U's for observations
        U1q       = Cov_ll_ii*Bm_i*W_gg_ii * U1qc * Bm_i'; %idempotent, rank |r_U|
        U1mq(i,:) = U1q(:)';
        U2        = Cov_ll_ii*Bm_i*W_gg_ii * U2c * Bm_i'; %idempotent, rank U-|r_U|
        U2m(i,:)  = U2(:)';
        % check     = U1q + U2 + Rii;   % idempotent, rank G
        
        % sensitivity factors for seleced parameters
        muv1(i)   = sqrt(real(max(eig(U1q * Rii_inv))));
    end
end
