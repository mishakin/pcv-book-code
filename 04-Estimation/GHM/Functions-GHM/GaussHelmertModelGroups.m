%% Gauss Helmert Model for Groups
%
% g_i(^l_i,^x) = 0  group of constraints
%
% all observational groups have the same dimension d_I
% all observational groups are mutually uncorrelated
% all groups of constraints have the same size d_G
% constraints must not overlap, ie must not use the same observation
% covariance matrix of estimates is assumed to be full
% Only one vector of unknowns
% no robust estimation
%
% lm       = I x d_I matrix of observational groups
% Cov_ll_m = I x d_i^2 matrix of vectorized covariance matrices
% cg_f     = d_G x 1 constraint function -> cg, A, B
% xa       = Ux1 vector of approximate parameters
% ux       = function to update parameters x = ux(xa,dx)
% sx       = Ux1 vector of standard deviations of parameters
% Tx       = threshold for convergence
% maxiter  = maximum number of iterations
% rU       = range of unknown parameters of interest
%
% est_x    = Ux1 estimated parameter
% Cov_xx   = UxU theoretical covariance matrix of estimated parameters
%                may be full
% sigma_0q = estimated variance factor (=1 if R=0)
% R        = redundancy
% vv       = Ixd_I vector of estimated corrections
% Xv       = Ix1 vector of Xhi^2-test statistic (using sigma_0=1)
% Rim      = I x d_I^2 redundancy matrices as row vectors
% nabla_lv = Ix1 minimum detectabele outlier
% muv      = Ix1 vector of senstivity factors wrt to all parameters
% muv1     = Ix1 vector of senstivity factors wrt to selected group
% Uv1q     = I x d_I^2 matrices of contributions to selected group
% Uv2      = I x d_I^2 matrices of contribution to remaining parameters
%
% [est_x,Cov_xx,sigma_0q,R,vv,Xv,Rim,nabla_lv,muv,muv1,uv1q,uv2]=...
%     GaussHelmertModelGroups(lm,Cov_ll_m,cg_f,xa,ux,sx,Tx,maxiter,rU)
%
% Wolfgang Förstner 2016-06-06
% wfoerstn@uni-bonn.de

function [est_x,Cov_xx,sigma_0q,R,vv,Xv,Rim,nabla_lv,muv,muv1,Uv1q,Uv2]=...
    GaussHelmertModelGroups(lm,Cov_ll_m,cg_f,xa,ux,sx,Tx,maxiter,rU)

%% initialization ---------------------------------------------------------
% numbers N, G, and U
% numbers I, d_I
[I,d_I]    = size(lm);
[cg,Am,Bm] = cg_f(lm(1,:)',lm(1,:)',xa);
nz_A = size(Am,2);
nz_B = size(Bm,2);
d_G        = size(cg,1);
N          = size(Cov_ll_m,2);
G          = I * d_G;
U          = size(xa,1);

% redundancy
R = G-U;
if R < 0
    disp('not enough observations')
    return;
end
%% initialization ---------------------------------------------------------
nu  = 0;                    % index of iterations
la  = lm;                   % matrix of approximate fitted observations
lnu = la;                   % matrix of current fitted observations
xnu = xa;                   % vector of approxmate parameters
s   = 0;                    % iteration indicator

Am     = spalloc(G,U,G*d_G*nz_A);
Bmt    = spalloc(G,N,G*d_G*nz_B);
Cov_ll = spalloc(N,N,N*d_I^2);
W_gg   = spalloc(G,G,G*d_G^2);

for iter = 1: maxiter
    
    for i=1:I
        % residual constraints, Jacobains
        [cg_i,Amt_i,Bmt_i] = cg_f(lm(i,:)',lnu(i,:)',xnu);
        Am((i-1)*d_G+1:i*d_G,:)                 = Amt_i;                    %#ok<*SPRIX>
        Bmt((i-1)*d_G+1:i*d_G,(i-1)*d_I+1:i*d_I)= Bmt_i;
        cg((i-1)*d_G+1:i*d_G)                   = cg_i;
        Cov_ll_i                                = reshape(Cov_ll_m(i,:),d_I,d_I);
        % weight matrix of lv
        Cov_ll((i-1)*d_I+1:i*d_I,(i-1)*d_I+1:i*d_I) ...
            = Cov_ll_i;
        % weights of constraints
        W_gg_i     = inv(Bmt_i * Cov_ll_i * Bmt_i');
        W_gg((i-1)*d_G+1:i*d_G,(i-1)*d_G+1:i*d_G)  = W_gg_i;
    end
    Bm=Bmt';
    ABWB     = Am'  * W_gg;
    Nm       = ABWB * Am;               % normal equation matrix
    nv       = ABWB * cg;               % right hand sides
    
    % dNmis    = diag(1./sqrt(diag(Nm))); % condition Nm
    % T_Nm     = dNmis * Nm * dNmis;
    lambda_N = eigs(Nm);
    if abs(log(lambda_N(1)/lambda_N(U))) > 10
        disp('normal equation matrix nearly singular')
        return;
    end
    Cov_xx  = inv(Nm);                  % CovM of parameters
    delta_x = Nm \ nv;                  % correction of parameters
    
    nu = nu+1;
    
    % check for final iteration
    if max(abs(delta_x(:)./sx(:))) < Tx || nu == maxiter
        s = 2;
    end
    % correction of parameters
    xnu = ux(xnu,delta_x);
    
    % dx   = delta_x'
    % nu_x   = [nu,xnu']
    
    vv        = lnu - lm;                          % v^a
    vv_vector = reshape(vv',I*d_I,1);
    % residual constraints, Jacobains
    delta_l   = Cov_ll * Bm * W_gg ...
        * (cg - Am * delta_x) - vv_vector;
    % delta_ls  = delta_l';
    lnu       = lnu + reshape(delta_l,d_I,I)';     % fitted observations
    % correction to fitted observations
    vv        = lnu - lm;                          % residual matrix
    vv_vector = reshape(vv',I*d_I,1);
    % Omega
    % from residuals of constraints,
    % also works if covariance matrix of observations is singular
    cg = cg - Am * delta_x;
    Omega = cg' * W_gg * cg;
    
    if s==2
        break;
    end
end
if R > 0
    sigma_0q = Omega/R;
else
    sigma_0q = 1;
end
est_x = xnu;


%% diagnostics (useful if observations are groupwise uncorrelated) --------

[~,Xv,Rim,nabla_lv,muv,muv1,Uv1q,Uv2] = ...
    diagnostics_GHM_constraints_multi_d...
    (rU,I,d_I,d_G,Am,Bm,Cov_xx,Nm,W_gg,vv_vector);
