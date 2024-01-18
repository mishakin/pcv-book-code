%% diagnostics GMM 1d
%
% r_U      = range of paraemters of interest
% Am       = design matrix
% Cov_ll   = CovM of observations
% Cov_xx   = CovM of parameters
% W_xx     = weight matrix of parameters
% vv       = vector of residuals
%
% riv      = redundancy numbers of observations
% zv       = standardized resiuals
% nabla_lv = minimum size of detectable outliers
% muv      = sensitivity factors wrt to all parameters
% uv1q     = effect onty parameters of interest
% uv2      = effect onto nuisance parameters
%
% Wolfgang Förstner 6/6/2015
% wfoerstn@uni-bonn.de

function [riv,zv,nabla_lv,muv,muv1,uv1q,uv2] = ...
    diagnostics_GHM_1d(r_U,Am,Bm,Cov_ll,Cov_xx,W_xx,W_gg,vv)

delta_0 = 4.13;

U  = size(Cov_xx,1);

%% observations: detectatbility, statistical test
Rm     = Cov_ll*Bm*W_gg*Bm'-Cov_ll*Bm*W_gg*Am*Cov_xx*Am'*W_gg*Bm'; % redundancy matrix
riv    = diag(Rm);                         % redundancy numbers
zv     = vv./sqrt(riv.*diag(Cov_ll)+eps);  % normalized residuals
nabla_lv = delta_0 * sqrt(1./riv);

%% sensitivity wrt all parameters
muv    = sqrt((1-riv)./(riv+eps));         % sensitivity factors

%% sensitivity wrt selected parameters
if ~isempty(r_U)
    % ranges
    r_C= r_U;
    r_D= setdiff(1:U,r_U);
    % partitioned Jacobian matrix A=[C,D], G x [U_C,U_D]
    Cm = Am(:,r_C);         % G * U_C
    Dm = Am(:,r_D);         % G * U_D
    % reduced design matrix
    Cmr = Cm - Dm * inv(W_xx(r_D,r_D)) * W_xx(r_D,r_C);  %#ok<MINV> % C_reduced
    Cov_11 = Cov_xx(r_C,r_C);
    % effect of constraints onto parameters
    U1qc     = Cmr * Cov_11        * Cmr' * W_gg;         % U_C_bar
    U2c      = Dm  * inv(W_xx(r_D,r_D)) * Dm'  * W_gg;    %#ok<MINV> % U_D
    
    % ---  for checking
    % uv1qc    = diag(U1qc);                                 % u_C_bar
    % uv2c     = diag(U2c) ;                                 % u_D
    % Rmc      = eye(G) - Am*Cov_xx*Am'*W_gg;
    % rivc     = diag(Rmc);
    % check =uv1qc+uv2c+rivc ; % should be ones
    
    % sensitivity factors for constraints
    % wrt observations
    U1q     = Cov_ll*Bm*W_gg * U1qc *Bm';
    uv1q    = diag(U1q);
    U2      = Cov_ll*Bm*W_gg * U2c * Bm';
    uv2     = diag(U2);
    
    % ---  for checking
    % muv1c    = sqrt(uv1qc ./ (rivc+eps));
    % check =uv1q+uv2+riv; % should be ones
    
    % sensitivity factors for observations
    muv1    = sqrt(uv1q ./ (riv+eps));
end
