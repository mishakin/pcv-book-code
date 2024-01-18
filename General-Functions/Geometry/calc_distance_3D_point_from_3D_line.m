%% distance 3D point from 3D line
%
% Usage:
%   d = calc_distance_3D_point_from_3D_line(X,L)
%
%   X  - 4 x 1 homogeneous vector of 3D point
%   L  - 6 x 1 Plücker coordinates 
%
%   d  - Euclidean distance, inf if distance infinite
%
% Wolfgang Förstner 2010/02/14
% wfoerstn@uni-bonn.de 
%
% See also calc_distance_3D_point_from_linesegment, calc_distance_3D_point_from_triangle

function d = calc_distance_3D_point_from_3D_line(X,L)

% homogeneous parts
Xh = X(4);
Lh = L(1:3);
% check for elements at infinity
if Xh == 0 || norm(Lh) == 0
    d = inf;
else
    % distance
    d = norm(calc_S(X(1:3))*Lh-Xh*L(4:6))/norm(Xh*Lh);
end