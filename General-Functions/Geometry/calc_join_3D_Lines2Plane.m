%% plane through two 3D-lines L and M 
% precondition: lines are coplanar and not equal
%
% Usage:
%   A = calc_join_3D_Lines2Plane(L,M)
%
%   L, M - 6x1 Pl�cker Lines
%   A    - 4x1 plane
%
% Wolfgang F�rstner
% wfoerstn@uni-bonn.de 
%
% See also calc_intersect_3D_Lines2Point, calc_intersect_3D_Lines2Point_finite
% calc_join_3D_Lines2Plane_finite

function A = calc_join_3D_Lines2Plane(L,M)

% matrix for constraints Gamma(L)*A=0, ...
T = [calc_Gamma(L);...
     calc_Gamma(M)];

% nullspace of T yields the target point X
[~,~,V] = svd(T);
A = V(:,4);