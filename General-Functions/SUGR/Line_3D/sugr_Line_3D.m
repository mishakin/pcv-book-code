%% Create 3D line
%
% Line_3D = sugr_Line_3D(a1,a2)
%
% Usage:
%    Line_3D = sugr_Line_3D(L)
%       6-vector, CovM = 0; needs not be a Plücker vector
%    Line_3D = sugr_Line_3D(L,Chh)
%       homgeneous 6-vector, needs not be a Plücker vector
%
% Line_3D = structure
% *      .h     = spherically normalized Plücker vector
% *      .Crr   = reduced covariance matrix
% *      .Jr    = null space of [.h,Dual(.h)]'
% *      .type = 4
%
% Bernhard Wrobel, Wolfgang Förstner 7/2012
% wfoerstn@uni-bonn.de
%
% See also sugr_minimal_3D_Line, sugr_construct_join_Line_3D, sugr_constrain_3D_Line
% sugr_Point_3D , sugr_Plane_3D


function Line_3D = sugr_Line_3D(a1,a2)

% default
Chh  = zeros(6);

switch nargin
    case 0
        
    case 1 % 1 argument
        
        % (1) Plücker vector, CovM == 0
        % spherically normalized, Plücker constraint
        Line_3D = sugr_minimal_3D_Line(a1,Chh);
        
    case 2   % 2 argument
        
        % (2) homogeneous vector
        Line_3D = sugr_minimal_3D_Line(a1,a2);
        
end

Line_3D.type=4;
