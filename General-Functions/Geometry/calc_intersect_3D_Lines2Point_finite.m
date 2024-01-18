%% intersection point X of two 3D-lines L and M 
% precondition: L and M are not parallel, otherwise X = [0 0 0 0]'
%
% see PCV (13.306)
%
% Usage:
%   X = calc_intersect_3D_Lines2Point_finite(L,M)
%
%   L, M - 6x1 Plücker Lines
%   X    - 4x1 homogeneous 3D point closest to L and M
%
% Wolfgang Förstner 7/2017
% wfoerstn@uni-bonn.de 
%
% See also calc_intersect_3D_Lines2Point, calc_join_3D_Lines2Plane, 
% calc_join_3D_Lines2Plane_finite

function X = calc_intersect_3D_Lines2Point_finite(L,M)

% take safest point at finity
GL = calc_Gamma_reduced(L);
GM = calc_Gamma_reduced(M);
[~,il] = max(abs(GL(:,4)));
[~,ml] = max(abs(GM(:,4)));
P = GL(il,1:3)'/GL(il,4);
Q = GM(ml,1:3)'/GM(ml,4);

% directions
R = L(1:3);
S = M(1:3);

% solving linear equation system
RS = [R,-S];
Nm = RS'*RS;

if abs(det(Nm)) > 10^(-10)
    % R not parallel to S: finite point
    lm = (Nm)\RS'*(Q-P);
    % take midpoint
    X = [(P+lm(1)*R+Q+lm(2)*S)/2;1];
else
    % point indefinite
    X = zeros(4,1);
end
