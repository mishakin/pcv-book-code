%% condition points, determine conditioning matrix
% see PCV (6.137), but sigma*sqrt(dim) instead of max
%
% Usage:
%   [xc,M] = condition_Points(x)
%  
%	x  - Nxd matrix of Cartesian coordinates
% 
%   xc - Nxd matrix of Cartesian coordinates of conditioned points
%   M  - d+1 x d+1 matrix of conditioning
%
% Wolfgang Förstner 2/2013
% wfoerstn@uni-bonn.de 
%
% See also uncondition_Points, condition_Homography, uncondition_Homography

function [xc,M] = condition_Points(x)

% number and dimension
[N,d] = size(x);

% mean and scale
xm = mean(x);
C  = cov(x);
sigma = sqrt(d)*sqrt(trace(C)/d);

% conditioning matrix
M = [eye(d)/sigma  -xm/sigma;...
     zeros(1,d)    1];
 
% condition points 
xch =  (M * [x ones(N,1)]')';

% Cartesian coordinates
xc = xch(:,1:d);


