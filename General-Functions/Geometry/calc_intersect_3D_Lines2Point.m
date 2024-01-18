%% intersection point X of two coplanar 3D-lines L and M 
% precondition: lines are not equal
%
% Usage:
%   X = calc_intersect_3D_Lines2Point(L,M)
%
%   L, M - 6x1 Plücker Lines
%   X    - 4x1 homogeneous 3D point
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% See also calc_intersect_3D_Lines2Point_finite
% calc_join_3D_Lines2Plane, calc_join_3D_Lines2Plane_finite

function X = calc_intersect_3D_Lines2Point(L,M)

% matrix for constraints Gammadual(L)*X=0, ...
T = [calc_Gammadual(L);...
     calc_Gammadual(M)];

% nullspace of T yields the target point X
[~,~,V] = svd(T);

X = V(:,4);