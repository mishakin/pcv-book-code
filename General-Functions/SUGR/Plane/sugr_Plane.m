%% Create Plane
%
% Usage:
%
%    Plane = sugr_Plane(Ah)
%        Ah = 4x1 homogeneous coordinate vector, CovM = 0
%
%    Plane = sugr_Plane(n, d)
%        Hessian parameters, CovM = 0
%        n = 3x1 normal
%        d = 1x1 distance
%
%    Plane = sugr_Plane(Ah,Ahh)
%        Ah = 4x1 homogeneous coordinate vector
%        Ahh = 4x4 cov.matrx
%
%
% Plane = structure
% *           .h    = spherically normalized homogeneous coordinates
% *           .Crr  = reduced covariance matrix of l.h
% *           .type = 5
%
%
% Wolfgang Förstner 7/2012
% wfoerstn@uni-bonn.de
%
% See also sugr_Point_3D, sugr_Line_3D

function Plane = sugr_Plane(a1,a2)

% default
Chh  = zeros(4);


switch nargin
    case 1 % 4-vector
        %% (1) homogeneous vector, CovM == 0
        % minimal representation
        Plane   = sugr_minimal_vector(a1,Chh);
        
    case 2
        size_a1=size(a1,1);
        switch size_a1
            case 3
                %% (2) normal, distance, CovM == 0
                % minimal representation
                Plane   = sugr_minimal_vector([a1;-a2],Chh);
                
            case 4
                %% (2) homogeneous coordinate vector
                % minimal representation
                Plane   = sugr_minimal_vector(a1,a2);
        end
        
end

Plane.type = 5;
