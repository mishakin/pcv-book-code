% Abstand Punkt T und Dreieck D(X,Y,Z)
%
% Usage:
%   dist = calc_distance_3D_point_from_triangle(T,X,Y,Z)
%
%   T, X, Y, Z - 3D points, euclidean
%   dist       - double, shortest distance from T to triangle X-Y-Z 
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% last changes
% Susanne Wenzel 12/17
% wenzel@igg.uni-bonn.de
%
% See also calc_distance_3D_point_from_3D_line, calc_distance_3D_point_from_linesegment

function dist = calc_distance_3D_point_from_triangle(T,X,Y,Z)

% normal of triangle
V  = cross(X-Y,Z-Y);
Vu = [V;0];

% check geometric relation of given points
f = det([[X;1] [Y;1] [Z;1] Vu])+eps;
u = det([[Y;1] [Z;1] [T;1] Vu])/f;
v = det([[Z;1] [X;1] [T;1] Vu])/f;
w = det([[X;1] [Y;1] [T;1] Vu])/f;

if u >= 0 && v >= 0 && w >= 0
    % foot point of T lies inside X-Y-Z
    % distance point plane
    Am = [[X;1] [Y;1] [Z;1]];
    A = null(Am');
    A = A(1:4,1);
    dist = abs(T'*A(1:3)+A(4))/norm(A(1:3));   % see PCV (7.99), given T_h = 1
    % following PCV (7.99) it should look like
    % dist = abs([T;1]'*A)/norm(1*A(1:3));
end
if u <= 0 && v >= 0 && w >= 0
    % foot point of T lies outside triangel X-Y-Z, but next to 
    % line segment Y-Z
    dist = calc_distance_3D_point_from_linesegment(T,Y,Z);
end
if u >= 0 && v <= 0 && w >= 0
    % foot point of T lies outside triangel X-Y-Z, but next to 
    % line segment X-Z
    dist = distanz_Punkt_Kante(T,X,Z);
end
if u >= 0 && v >= 0 && w <= 0
    % foot point of T lies outside triangel X-Y-Z, but next to 
    % line segment X-Y
    dist = calc_distance_3D_point_from_linesegment(T,X,Y);
end
if u >= 0 && v <= 0 && w <= 0
    % foot point of T lies outside triangel X-Y-Z, 
    % outside linesegments
    % but next to point X
    dist = norm(T-X);
end
if u <= 0 && v >= 0 && w <= 0
    % foot point of T lies outside triangel X-Y-Z, 
    % outside linesegments
    % but next to point Y
    dist = norm(T-Y);
end
if u <= 0 && v <= 0 && w >= 0
    % foot point of T lies outside triangel X-Y-Z, 
    % outside linesegments
    % but next to point Z
    dist = norm(T-Z);
end

