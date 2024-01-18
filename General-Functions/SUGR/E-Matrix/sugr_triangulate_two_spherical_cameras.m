%% sugr_triangulate_two_spherical_cameras 
%
% Alg. 21: Optimal triangulation
%
% [X,sigma0,est_u,est_v, f] = sugr_triangulate_two_spherical_cameras(...
%                   b,R,u,v,sigma,Tol, max_iter, k, print_option)
%
% input:
% b,Rot    = relative orientation
% u,v      = 3x1 directions, spherically normalized
% sigma    = radial error of direction
% Tol      = tolerance for stopping iterations
% max_iter = maximum iterations
% k        = critical value
% print_option = boolean, deciding on printout
%
% output:
% X       = 4x1 vector, 3d point
% sigma0  = estimated sigma_0
% est_u, est_v = fitted observations
% f            = 0 at infinity
%              = 1 finite
%              = 2 backward
%              = 3 divergent
%
% Wolfgang Förstner 2014/04/09
% wfoerstn@uni-bonn.de 

function [X_est, est_sigma_0, est_u, est_v, f] = ...
    sugr_triangulate_two_spherical_cameras(...
    b, Rot, u, v, sigma, Tol, max_iter, k, print_option)
      
% essentiaol matrix
E = calc_S(b)*Rot';

% critical value for constraint
sigma_c= sigma*sqrt(u'*(E*E')*u+v'*(E'*E)*v);
crit_E = k*sigma_c;

cgl = u'*E*v;
if print_option > 0
    cgl=cgl                                                                %#ok<ASGSL,NOPRT>
end
% check for coplanarity
if abs(cgl) > crit_E
    X_est.h=zeros(4,1);
    X_est.Crr=zeros(3);
    X_est.type=3;
    est_sigma_0 = 0;
    f=2;
    est_u=u;
    est_v=v;
    return
end

% initialize
f  = 0;
ua = u;
va = v;
est_sigma_0=0;

% iterate
for nu = 1:max_iter
    g = ua'*E*va;               % constraint
    if print_option > 0
        nu=nu                                                              %#ok<FXSET,ASGSL,NOPRT>
        g=g                                                                %#ok<ASGSL,NOPRT>
    end
    % check for convergence
    if abs(g) < Tol*sigma
        break
    end
    J1 = null(ua');              % reducing Jacobians
    J2 = null(va');
    l  = [J1'*u;J2'*v];          % reduced observations
    B  = [J1'*E*va;J2'*E'*ua];   % Jacobian of constraint
    cg = -cgl;                   % residual constraint
    n = B'*B;
    Deltal = l + cg/n*B;         % corrections to reduced observations
    ua = ua+J1*Deltal(1:2);      % updated approximate values
    va = va+J2*Deltal(3:4);
    ua = ua/norm(ua);
    va = va/norm(va);
    est_sigma_0 = abs(cgl)/sigma_c; %/sqrt(n);  % estimated sigma_0
end
% determine homgeneous coordinates
est_u = ua;                     % final estimates
est_v = va;
est_w = Rot'*est_v;               % v in left system
if print_option > 0
    uvw= [est_u', est_v', est_w']                                          %#ok<NASGU,NOPRT>
end
m = calc_S(calc_S(b)*est_u)*b;
m = m /norm(m);                 % auxiliary vector
D = det([b, m, calc_S(est_u)*est_w]);
rs = m' * [est_w, est_u];
X = [rs(1)* est_u; D];
%  Xe   = X'
%  uvw  = [est_u,est_v,est_w]
%  mDrs = [m',D,rs]
X_est = sugr_cov_matr_X(X,sigma,est_u,est_w,rs);

% check for infinity


% critical value for Determinant
sigma_D = 2*sigma;         % conservative approximation
crit_D  = k*sigma_D;
if abs(D) < eps
    if print_option > 0
        disp('point at infinity')
    end    
    return
end
d=rs/D;
% check for backward direction and infinity
if sign(d(1)) < 0 || sign(d(2)) < 0
    if abs(D) < crit_D
        % point at infinity
        if print_option > 0
            disp('point at infinity')
        end
        
        return
    else
        % rays divergent
        f=1;
        X_est.h=zeros(4,1);
        X_est.Crr=zeros(3);
        X_est.type=3;
        if print_option > 0
            disp('backward point')
        end
        %dbstop
        return
    end
else
    if print_option > 0
        disp('signs ok')
    end
end



