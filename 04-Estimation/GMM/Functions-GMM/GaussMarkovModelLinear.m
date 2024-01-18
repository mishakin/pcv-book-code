%% Gauss-Markov model linear, individual observations
%
% lv       = Nx1 vector of observations
% Cov_ll   = NxN covariance matrix of observations (Cov_ll^a)
% Am       = NxU Jacobian
% av       = Nx1 additive vector
% rU       = range of unknown parameters of interest
%
% est_x    = Ux1 estimated parameter
% Cov_xx   = UxU theoretical covariance matrix of estimated parameters
% sigma_0q = estimated variance factor (=1 if R=1)
% R        = redundancy
% vv       = Nx1 vector of estimated corrections
% zv       = Nx1 vector of standardized residuals (using sigma_0=1)
% riv      = Nx1 vector of redundacy numbers
% mu       = Nx1 vector of senstivity factors
%
% [est_x,Cov_xx,sigma_0q,R,vv,zv,riv,muv]=...
%            GaussMarkovModelLinear(l,Cov_ll,A,a)
%
% Wolfgang Förstner 2015-06-04
% wfoerstn@uni-bonn.de


function [est_x,Cov_xx,sigma_0q,R,vv,zv,riv,nabla_lv,muv,muv1,uv1q,uv2]=...
    GaussMarkovModelLinear(lv,Cov_ll,Am,av,rU)

%% initialization

% number N and U
[N,U] = size(Am);

% redundancy
R = N-U;
if R < 0
    disp('not enough observations')
    return;
end

%% estimation
W_ll   = inv(Cov_ll);                  % weight matrix of lv
Bm     = W_ll*Am;                         %#ok<*MINV> % ancillary matrix
Nm     = Bm'*Am;                       % normal equation matrix
mv     = Bm'*(lv-av);                  % right hand sides
Cov_xx =inv(Nm);                       % covariance matrix of est. parameters
est_x  = Cov_xx*mv;                    % estimated parameters
vv     = Am*est_x+av-lv;               % estimated residuals
if R > 0
    sigma_0q = vv'*W_ll*vv/R;          % estimated variance factor
else
    sigma_0q = 1;
end

%% diagnostics (useful if observations are uncorrelated)
[riv,zv,nabla_lv,muv,muv1,uv1q,uv2] = diagnostics_GMM_1d(rU,Am,Cov_ll,W_ll,Cov_xx,Nm,vv);
