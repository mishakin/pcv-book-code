%% plane parallel to and through midpoint between two 3D-lines L and M 
% precondition: lines not parallel, otherwise A = [0 0 0 0]'
%
% Usage:
%   A = calc_join_3D_Lines2Plane_finite(L,M)
%
%   L, M - 6x1 Plücker Lines
%   A    - 4x1 plane
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% See also calc_intersect_3D_Lines2Point, calc_intersect_3D_Lines2Point_finite
% calc_join_3D_Lines2Plane_finite

function A = calc_join_3D_Lines2Plane_finite(L,M)

% determine point closest to L and M
X = calc_intersect_3D_Lines2Point_finite(L,M);

% normal on L and M
N = cross(L(1:3),M(1:3));

if norm(X) > 10^(-10) && norm(N) > 10^(-10) 
    % finite plane
    A = [N*X(4); -N'*X(1:3)];
else
    % indefinite plane
    A = zeros(4,1);
end