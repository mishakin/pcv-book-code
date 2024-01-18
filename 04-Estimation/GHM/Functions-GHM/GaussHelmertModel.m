%% Gauss Helmert Model
%
% g(^l,^x) = 0
%
% lv       = Nx1 vector of observations
% Cov_ll   = NxN covariance matrix of observations (Cov_ll^a)
% cg_f     = Gx1 constraint function -> cg, A, B
% xa       = Ux1 vector of approximate parameters
% sx       = Ux1 vector of standard deviations of parameters
% Tx       = thershold for convergence
% maxiter  = maximum number of iterations
% rU       = range of unknown parameters of interest
%
% est_x    = Ux1 estimated parameter
% Cov_xx   = UxU theoretical covariance matrix of estimated parameters
%                may be full
% sigma_0q = estimated variance factor (=1 if R=0)
% R        = redundancy
% vv       = Nx1 vector of estimated corrections
% zv       = Nx1 vector of standardized residuals (using sigma_0=1)
% riv      = Nx1 vector of redundancy numbers
% nabla_lv = minimum detectabel outlier
% muv      = Nx1 vector of senstivity factors wrt to all parameters
% muv1     = Nx1 vector of senstivity factors wrt to selected group
% uv1q     = Nx1 vector of contributions to selected group
% uv2      = Nx1 vector of contribution to remaining parameters
%
% [est_x,Cov_xx,sigma_0q,R,vv,zv,riv,nabla_lv,muv,muv1,uv1q,uv2]=...
%   GaussHelmertModelNonlinearNotrobust(lv,Cov_ll,cg_f,xa,sx,Tx,maxiter,rU)
%
% Wolfgang Förstner 2016-06-06
% wfoerstn@uni-bonn.de


function [est_x,Cov_xx,sigma_0q,R,vv,zv,riv,nabla_lv,muv,muv1,uv1q,uv2] = ...
    GaussHelmertModel(lv,Cov_ll,cg_f,xa,sx,Tx,maxiter,rU)

%% initialization of estimation

% numbers G, and U
[cg,Am,Bm] = cg_f(lv,lv,xa);
G = size(cg,1);
U = size(xa,1);

% redundancy
R = G-U;
if R < 0
    disp('not enough observations')
    return;
end
%% initialization of iterations

nu  = 0;                         % number of iterations
la  = lv;
lnu = la;
xnu = xa;
s   = 0;

%% iterations

for iter = 1: maxiter
    [cg,Am,Bm] = cg_f(lv,lnu,xnu);      % residual if constraints, Jacobians
    %     W_ll     = inv(Cov_ll);             % weight matrix of lv
    W_gg     = inv(Bm' * Cov_ll * Bm);
    ABWB     = Am' * W_gg;                                                  %#ok<MINV>
    Nm       = ABWB * Am;               % normal equation matrix
    nv       = ABWB * cg;               % right hand sides
    lambda_N = eigs(Nm);
    if abs(log(lambda_N(1)/lambda_N(U))) > 10
        disp('normal equation matrix nearly singular')
        return;
    end
    Cov_xx = inv(Nm);                      % CovM of parameters
    delta_x = Cov_xx * nv;                 %#ok<MINV>  % correction of parameters
    
    nu = nu+1;
    
    % check for final iteration
    if max(abs(delta_x(:)./sx(:))) < Tx || nu == maxiter
        s = 2;
    end
    % correction of parameters
    xnu = xnu + delta_x;
    
    disp([num2str(nu),'. Iteration: delta_x = [',num2str(delta_x'),...
        '], est_x = [',num2str(xnu'),']'])
    
    % correction to fitted observations
    vv      = lnu - lv;             % residuals
    delta_l = Cov_ll * Bm * W_gg * (cg - Am * delta_x) - vv;                %#ok<MINV>
    lnu     = lnu + delta_l;        % fitted observations
    
    if s==2
        break;
    end
end
if R > 0
    sigma_0q = (cg' * W_gg *cg)/R;                                          %#ok<MINV>
else
    sigma_0q = 1;
end
est_x = xnu;

%% diagnostics (useful if observations are uncorrelated)
[riv,zv,nabla_lv,muv,muv1,uv1q,uv2] = ...
    diagnostics_GHM_1d(rU,Am,Bm,Cov_ll,Cov_xx,Nm,W_gg,vv);
