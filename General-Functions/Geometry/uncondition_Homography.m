%% uncondition Homography
% see PCV (6.138)
%
% Usage:
%   H = uncondition_Homography(Hc, Tl, Tr)
%  
%   Hc  - homogeneous conditioned Homography
%   Tl - condition matrix point set 1
%   Tr - condition matrix point set 2
% 
%   H - Homography, unconditioned
%
% Wolfgang Förstner 2/2013
% wfoerstn@uni-bonn.de 
%
% See also condition_Points, uncondition_Points, condition_Homography

function H = uncondition_Homography(Hc, Tl, Tr)

% unconditioned Homography
H =  inv(Tl) * Hc * Tr;                                                     %#ok<MINV>


