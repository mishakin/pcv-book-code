%% Create 3D point
%
% Point_3D = sugr_Point_3D(a1,a2,a3,a4,a5,a6)
%
% Usage:
% (1)   Point_3D = sugr_Point_3D(Xe)
%         Euclidean vector, CovM = 0
%
% (1)   Point_3D = sugr_Point_3D(Xh)
%         homogeneous vector, CovM = 0
%
% (2)   Point_3D = sugr_Point_3D(Xe,Cee)
%         Euclidean vector and CovM
%
% (2)   Point_3D = sugr_Point_3D(Xh,Chh)
%         homgeneous vector and CovM
%
% (6)   Point_3D = sugr_Point_3D(X,Y,Z,sX,sY,sZ)
%         independent Euclidean coordinates with standard dev
%
% * Point_3D = structure
%      .h     = spherically normalized homogeneous coordinates
%      .Crr   = reduced covariance matrix
%      .Jr    = null space of .h'
%      .type  = 3
%
% Bernhard Wrobel, Wolfgang Förstner 2/2013
% wfoerstn@uni-bonn.de
%
% See also sugr_Line_3D, sugr_Plane_3D, sugr_Point_2D


function Point_3D = sugr_Point_3D(a1,a2,a3,a4,a5,a6)

% default
Chh  = zeros(4);

switch nargin
    case 0
    case 1 % 3- or 4-vector
        size_a1 = size(a1,1);
        switch size_a1
            case 3
                %% (1) Euclidean vector, CovM == 0
                % spherically normalized
                Point_3D = sugr_minimal_vector([a1;1],Chh);
                
            case 4
                %% (1) homogeneous vector, CovM == 0
                Point_3D = sugr_minimal_vector(a1,Chh);
                
        end
    case 2
        
        size_a2 = size(a1,1);
        switch size_a2
            case 1
                Point_3D = sugr_minimal_vector([a1(:);1],[a2,zeros(3,1);zeros(1,4)]);
                
            case 3
                %% (2) Euclidean coordinates, CovM == 0
                Point_3D = sugr_minimal_vector([a1;1],[a2,zeros(3,1);zeros(1,4)]);
                
            case 4
                %% (2) homogeneous vector
                % spherically normalized
                Point_3D = sugr_minimal_vector(a1,a2);
        end
    case 6
        %% (6) homogeneous coordinates
        % individual standard deviations
        Point_3D = sugr_minimal_vector(...
            [a1;a2;a3;1],...
            [diag([a4,a5,a6]).^2,zeros(3,1);[0,0,0,1]]);
        
        
end

Point_3D.Jr = null(Point_3D.h');

Point_3D.type=3;

