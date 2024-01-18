%% Create 2D point
%
% Usage:
%    Point_2D = sugr_Point_2D(xe)
%         Euclidean vector, CovM = 0
%
%    Point_2D = sugr_Point_2D(xh)
%         homogeneous vector, CovM = 0
%
%    Point_2D = sugr_Point_2D(x,y)
%         Euclidean coordinates, CovM = 0
%
%   Point_2D = sugr_Point_2D(xe,Cee)
%         Euclidean vector with Cov.matrix
%
%    Point_2D = sugr_Point_2D(xh,Chh)
%          homgeneous vector with Cov.matrix
%
%    Point_2D = sugr_Point_2D(u,v,w)
%         homogeneous coordinates, CovM = 0
%
%    Point_2D = sugr_Point_2D(x,y,Cee)
%         Euclidean coordinates with Euclidean Cov.matrix
%
%    Point_2D = sugr_Point_2D(x,y,sx,sy)
%         independent Euclidean coordinates with standard deviations
%
%    Point_2D = sugr_Point_2D(u,v,w,Chh)
%         homogeneous coordinates with homogeneous Cov.matrix
% 
%  * Point_2D = structure 
%        .h    = spherically normalized homogeneous coordinates
%        .Crr  = reducedcovariance matrix
%        .Jr   = null space of .h'
%        .type = 1
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de
%
% See also sugr_Line_2D, sugr_Line_3D, sugr_Plane_3D, sugr_Point_3D

function Point_2D = sugr_Point_2D(a1,a2,a3,a4)

% default
Chh  = zeros(3);

switch nargin
    case 0
    case 1 % 2- or 3-vector
        size_a1 = size(a1,1);
        switch size_a1
            case 2 
%% (1) Euclidean vector, CovM == 0
                % spherically normalized
                Point_2D = sugr_minimal_vector([a1;1],Chh);
                
            case 3 
%% (1) homogeneous vector, CovM == 0
                Point_2D = sugr_minimal_vector(a1,Chh);
                
        end
    case 2

        size_a2 = size(a1,1);
        switch size_a2
            case 1 
%% (2) Euclidean coordinates, CovM == 0
                Point_2D = sugr_minimal_vector([a1;a2;1],Chh);
                
            case 2
%% (2) Euclidean vector
                % spehrically normalized 
                Point_2D = sugr_minimal_vector([a1;1],[a2 zeros(2,1); zeros(1,3)]);

            case 3
%% (3) uncertain homogeneous coordinates
                Point_2D   = sugr_minimal_vector(a1,a2);
        end
    case 3
        size_a3 = size(a3,1);
        switch size_a3
            case 1 
%% (3) homogeneous vector, CovM == 0
                 % spehrically normalized 
                Point_2D = sugr_minimal_vector([a1;a2;a3],Chh);
                
            case 2
%% (3) Euclidean coordinates
                % spehrically normalized 
                Point_2D = sugr_minimal_vector([a1;a2;1],[a3 zeros(2,1); zeros(1,3)]);
        end
        
    case 4
        size_a4 = size(a4,1);
        switch size_a4
            case 1 
%% (4) Euclidean vector, independent coordinates
                % spehrically normalized 
                Point_2D = sugr_minimal_vector([a1;a2;1],...
                           [ a3^2, 0, 0; 0, a4^2, 0; zeros(1,3)]);
                
            case 3
%% (4) homogeneous coordinates  
                % spehrically normalized 
                Point_2D = sugr_minimal_vector([a1;a2;a3],a4);
        end
        
end

Point_2D.Jr = null(Point_2D.h');
Point_2D.type=1;

