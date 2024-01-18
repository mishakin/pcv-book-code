%% Create pose
%
% Usage:
%
%    Point_Projection_3D_2D = sugr_Projection_3D_2D(P)
%        3x4 matrix, CovM = 0
%
%    Point_Projection_3D_2D = sugr_Projection_3D_2D(P,CovM)
%        3x4 matrix, 12x12 covariance matrix of vec(P)
%
%    Point_Projection_3D_2D = sugr_Projection_3D_2D(K,R,Z)
%        3x3-matrices, 3-vector, CovM = 0
%
% Point_Projection_3D_2D = structure
% *             .P   = 3 x 4-Matrix, Frobenius norm = 1, det P(1:3,1:3) > 0
% *             .Crr = vec of covariance matrix of 11 reduced parameters
% *             .type = 20
%
% Wolfgang Förstner 2/2013
% wfoerstn@uni-bonn.de
%
% See also sugr_minimal_Point_Projection_3D_2D, sugr_minimal_vector

function Point_Projection_3D_2D = sugr_Projection_3D_2D(a1,a2,a3)

% default
Chh  = zeros(12);
%
switch nargin
    
    case 1
        % 3x4-matrix
        P = a1*sign(det(a1(1:3,1:3)));
        Point_Projection_3D_2D = ...
            sugr_minimal_Point_Projection_3D_2D(P,Chh);
        
    case 2
        % P, Cpp
        P = a1*sign(det(a1(1:3,1:3)));
        Point_Projection_3D_2D = ...
            sugr_minimal_Point_Projection_3D_2D(P,a2);
        
    case 3 %
        % K, R, Z
        P = a1*a2*[eye(3),-a3];
        P = P*sign(det(P(1:3,1:3)));
        Point_Projection_3D_2D = ...
            sugr_minimal_Point_Projection_3D_2D(P,Chh);       
        
end

Point_Projection_3D_2D.type = 134;

