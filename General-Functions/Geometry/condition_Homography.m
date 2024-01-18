%% condition Homography
% see PCV (6.138)
%
% Usage:
%   Pc = condition_Homography(H, Tl, Tr)
%  
%   H  - homogeneous Homography
%   Tl - condition matrix point set 1
%   Tr - condition matrix point set 2
% 
%   Hc - Homography, conditioned
%
% Wolfgang Förstner 2/2013
% wfoerstn@uni-bonn.de 
%
% See also condition_Points, uncondition_Points, uncondition_Homography

function Hc = condition_Homography(H, Tl, Tr)

% conditioned Homography
Hc =  Tl * H * inv(Tr);                                                     %#ok<MINV>


