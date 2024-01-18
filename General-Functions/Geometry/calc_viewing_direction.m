%% Viewing direction from P of perspective camera
%
% Usage:
%   d = calc_viewing_direction(P)
%  
%   P - 4 x 3-projection matrix P = K*R*[I|-Z] 
%
%   d - 3x1 vector viewing direction, not normalized
%
% Wolfgang Förstner 4/2013
% wfoerstn@uni-bonn.de 
%
% See also calc_P_from_KRZ, calc_KRZ_from_P, calc_Q_from_P

function d = calc_viewing_direction(P)

d = -P(3,1:3)'*det(P(:,1:3));
