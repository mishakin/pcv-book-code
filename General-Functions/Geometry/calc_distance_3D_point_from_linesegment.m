%% distance of a point T to line segment K(X,Y)
%
% Usage:
%   dist = calc_distance_3D_point_from_linesegment(T,X,Y)
%
%   T, X, Y - 3D points, euclidean
%   dist    - double, shortest distance from T to segment X-Y 
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% See also calc_distance_3D_point_from_3D_line, calc_distance_3D_point_from_triangle

function dist = calc_distance_3D_point_from_linesegment(T,X,Y)

% direction
Lh = Y - X;
t = (T-X)'*Lh/norm(Lh)^2;
if t < 0
    dist = norm(T - X);
else
    if t>1
        dist = norm(T - Y);
    else
        L0   = cross(X,Y);
        dist = norm(cross( T,Lh)-L0)/norm(Lh);
    end
end

