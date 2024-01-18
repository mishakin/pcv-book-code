%% diagnostics 1d
%
% given all parameters of a linearized estimation problem
% determine the essential diagnostic values
% assuming the number N of observations is small enough
%    (requires several NxN matrices)
%
% r_U    = index set for parameters to be estimated (others are nuisance)
% Am     = Jacobian
% Cov_ll = CovM of observations
% W_ll   = WeightM of observations
% Cov_xx = CovM of all estimated parameters
% W_xx   = WeightM of all estimated parameters
% vv     = residuals
%
% riv    = redundancy numbers
% zv     = standardized test statistics
% nabla_lv = lowe bound for detectable outliers
% muv      = influence factor for all parameters
% muv1     = influence factor for parameters specified by r_U
% uv1q     = diagonal elements of \bar U_1 (numerator of muv1^2)
% uv2      = diagonal elements of U2
%
% Wolfgang Förstner 10/2016
% wfoerstn@uni-bonn.de

function [riv,zv,nabla_lv,muv,muv1,uv1q,uv2] = ...
    diagnostics_GMM_1d(r_U,Am,Cov_ll,W_ll,Cov_xx,W_xx,vv)

% non-centrality parameter for alpha_0 = 0.001, beta_0 = 0.8
delta_0 = 4.13;

[N,U]  = size(Am);
%% observations: detectatbility, statistical test
Rm     = eye(N)-Am*Cov_xx*Am'*W_ll;        % redundancy matrix
riv    = diag(Rm);                         % redundancy numbers
zv     = vv./sqrt(riv.*diag(Cov_ll)+eps);  % normalized residuals
nabla_lv = delta_0 * sqrt(1./riv);

%% sensitivity wrtall parameters
muv    = sqrt((1-riv)./(riv+eps));         % sensitivity factors

%% sensitivity wrt selected parameters
if ~isempty(r_U)
    % ranges
    r_C= r_U;
    r_D= setdiff(1:U,r_U);
    % partitioned design matrix A=[C,D]
    Cm = Am(:,r_C);
    Dm = Am(:,r_D);
    % reduced design matrix
    Cmr = Cm - Dm * inv(W_xx(r_D,r_D)) * W_xx(r_D,r_C);  %#ok<*MINV> % C_reduced
    Cov_11 = Cov_xx(r_C,r_C);
    % effect onto parameters
    U1q     = Cmr * Cov_11        * Cmr' * W_ll;         % U_C_bar
    uv1q    = diag(U1q);                                 % u_C_bar
    U2      = Dm  * inv(W_xx(r_D,r_D)) * Dm'  * W_ll;    % U_D
    uv2     = diag(U2) ;                                 % u_D
    % sensitivity factors
    muv1    = sqrt(uv1q ./ riv);
    
    %     check = uv1q+uv2+riv;  % should be ones
end
