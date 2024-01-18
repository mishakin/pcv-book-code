%% uncondition points, determine conditioning matrix
% see PCV (6.137)
%
% Usage:
%   x = uncondition_Points(xc,M)
%  
%   xc - Nxd matrix of Cartesian coordinates of conditioned points
%   M  - d+1 x d+1 matrix of conditioning
%
%	x  - Nxd matrix of Cartesian coordinates
%
% Wolfgang Förstner 2/2013
% wfoerstn@uni-bonn.de 
%
% See also condition_Points, condition_Homography, uncondition_Homography

function x = uncondition_Points(xc,M)

% number and dimension
[N,d]  =size(xc);

% condition points 
xh =  (inv(M) * [xc ones(N,1)]')';                                          %#ok<MINV>

% Cartesian coordinates
x = xh(:,1:d);


