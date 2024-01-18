%% determine CovM for triangulated point (cov_matr_X)
%
% determine CovM from estimated directions and 3D point 
% following PCV (13.268)
% directions are assumed to follwo epipolar constraint
%
% X_est = sugr_cov_matr_X(X,sigma,est_u,est_w,rs)
%
% X           = 4x1 homgeneous vector, spherically normaized
% sigma        = standard deviation of directions, isotropic
% est_u,       = estimated directions in first camera system
% est_w        = est_v in first camera system, w=R'*v
% rs           = 2x1 vector of relative distances
% 
% X_est        3D struct of point
% *     .h     = Xs
% *     .Crr   = CovM of reduced spherically homogeneous point coordinates
%                
% Wolfgang Förstner 
% wfoerstn@uni-bonn.de
%

function X_est = sugr_cov_matr_X(X,sigma,est_u,est_w,rs)

% determine basis
nv = calc_S(est_u)*est_w;
nv = nv/norm(nv);
rv = calc_S(nv)*est_u;
sv = calc_S(nv)*est_w;
C =  [inv(rv*rv'/rs(1)^2+sv*sv'/rs(2)^2+ ...
       nv*nv'*(rs(1)^2+rs(2)^2)/(rs(1)^2*rs(2)^2)+eye(3)*10^(-12)...
       )*sigma^2 ...
       zeros(3,1); [0,0,0,0]];                                              %#ok<MINV>
X_est = sugr_Point_3D(X,C);
