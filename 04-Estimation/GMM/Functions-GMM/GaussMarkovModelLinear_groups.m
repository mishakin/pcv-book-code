%% Gauss-Markov model linear, groups of uncorrel. observations
%
% all observational groups have the same dimension
% all observational groups are mutually uncorrelated
%
% input matrix Am may be sparse
% covariance matrix of estimates is assumed to be full
%
% lm       = I x d_I matrix of observational groups
% Cov_ll_m = I x d_i^2 matrix of vectorized covariance matrices
% Am       = NxU Jacobian (may be sparse)
% av       = Nx1 additive vector
% rU       = range of unknown parameters of interest
%
% est_x    = Ux1 estimated parameter
% Cov_xx   = UxU theoretical covariance matrix of estimated parameters
% sigma_0q = estimated variance factor (=1 if R=0)
% R        = redundancy
% vm       = I x d_I matrix of estimated residuals
% Xv       = Ix1 vector of test statistics  (Xhi^2_d_I using sigma_0=1)
% Rim      = I x d_i^2 matrix of redundacy matrices R_ii
% nabla_lv = size of minimal detectable errors
%            (from max eigenvalue of C_ll W_vv C_ll * delta_0
% muv      = Ix1 vector of senstivity factors
%               (mu_i:=mu_ix, cf. (4.307))
% muv1     = Ix1 vector of senstivity factors for selected parameters
%               (mu_i:=mu_ik, cf. (4.315))
% Um1q     = I x d_i^2 matrix of reduced U_ii-matrices
%              for selected parameters (cf. bar(U)_k(4.127))
% Um2      = I x d_i^2 matrix of U_ii-matrices
%              for not selected parameters (cf. U_t (4.127))
%
% [est_x,Cov_xx,sigma_0q,R,vm,Xv,Rim,muv]=...
%            GaussMarkovModelLinear_groups(lm,Cov_ll,Am,av)
%
% Wolfgang Förstner 2015-06-04
% wfoerstn@uni-bonn.de


function [est_x,Cov_xx,sigma_0q,R,vi,Xv,Rim,nabla_lv,muv,muv1,Um1q,Um2] = ...
    GaussMarkovModelLinear_groups(lm,Cov_ll_m,Am,av,rU)

%% initialization

% numbers I, d_I
[I,d_I]= size(lm);

% number N and U
[N,U] = size(Am);

% redundancy
R = N-U;
if R < 0
    disp('not enough observations')
    return;
end

% reshape observations and covariance matrix
lv     = zeros(N,1);
Cov_ll = spalloc(N,N,I*d_I^2);
for i = 1:I
    lv(d_I*i-(d_I-1):d_I*i)=lm(i,:)';
    Cov_ll(d_I*i-(d_I-1):d_I*i,d_I*i-(d_I-1):d_I*i) = ...
        reshape(Cov_ll_m(i,:),d_I,d_I);
end
%% estimation
W_ll   = inv(Cov_ll);                  % weight matrix of lv
Bm     = W_ll*Am;                      %#ok<*MINV> % ancillary matrix
Nm     = Bm'*Am;                       % normal equation matrix
mv     = Bm'*(lv-av);                  % right hand sides

% use sparseinv in Matlab, otherwise inv
%Cov_xx = sparseinv(Nm);               % covariance matrix of est. parameters
Cov_xx = inv(Nm);                      % covariance matrix of est. parameters
est_x  = Nm\mv;                        % estimated parameters
vv     = Am*est_x+av-lv;               % estimated residuals
if R > 0
    sigma_0q = vv'*W_ll*vv/R;          % estimated variance factor
else
    sigma_0q = 1;
end
%% diagnostics
[vi,Xv,Rim,nabla_lv,muv,muv1,Um1q,Um2]=...
    diagnostics_GMM_multi_d(rU,I,d_I,Am,Cov_ll,W_ll,Cov_xx,Nm,vv);
