%% robustly estimate sigma from residuals squared
%
% s = sugr_estimate_robust_scale_from_residuals(r,d)
%
% r = vector of absolute residuals
%     assumtption: at least 10 % inliers chi^2_d distributed
% d = degrees of freedom of residuals
% s = robust sigma estimate
%
% see adaptive least Kth-order squares estimator by Lee et al. (1998)
% adapted to magnitude of multi-dimensional residuals 
%
% Wolfgang Förstner 2/2012
% wfoerstn@uni-bonn.de 

function s0 = sugr_estimate_robust_scale_from_residuals(r,d)

N = length(r);                                % # data
p = 0.1;                                      % minimum percentage of inliers 

s0 = kth_element(r,N,p*N)/icdf('chi2',p,d);   % initial estimate

% b = 1.06*s0/N^(1/5);                          % bandwidth 


