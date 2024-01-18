%% create homogeneous 2D line
%
% Line_2D = sugr_Line_2D(a1,a2,a3,a4,a5)
%
% input cases
% (1) le             Hessian parameters [phi,d], CovM = 0
% (1) lh             homogeneous coordinate vector, CovM = 0
% (2) phi, d         Hessian parameters, CovM = 0
% (2) le,Cee         Euclidean coordinate vector
% (2) lh,Chh         homogeneous coordinate vector
% (3) a,b,c          homogeneous coordinates, CovM == 0
% (3) phi,d,Cee      Euclidean coordinates
% (4) u,v,w,Chh      homogeneous coordinates
% (5) x,y,phi,sp,sq  centroid form
%
% Line_2D = structure
% *           .h    = spherically normalized homogeneous coordinates
% *           .Crr  = reduced covariance matrix of l.h
% *           .Jr   = nullspace of .h'
% *           .type = 2
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de
%
% See also sugr_Line_2D_hom2Hes, sugr_Line_2D_hom2cen,
% sugr_Line_2D_Hes2hom, sugr_Line_2D_Hes2cen, sugr_Line_2D_cen2hom, 
% sugr_get_Euclidean_Line_2D, sugr_get_centroid_Line_2D

function Line_2D = sugr_Line_2D(a1,a2,a3,a4,a5)

% default
Chh  = zeros(3);

switch nargin
    case 1 % 2- or 3-vector
        size_a1=size(a1,1);
        switch size_a1
            case 2
                %% (1) Euclidean vector (angle-distance), CovM == 0
                % minimal representation
                Line_2D   = sugr_minimal_vector([cos(a1(1)); sin(a1(1)); -a1(2)],Chh);
                
            case 3
                %% (1) homogeneous vector, CovM == 0
                % minimal representation
                Line_2D   = sugr_minimal_vector(a1,Chh);
                
        end
        
    case 2
        size_a1=size(a1,1);
        switch size_a1
            case 1
                %% (2) Euclidean coordinates, CovM == 0
                % minimal representation
                Line_2D   = sugr_minimal_vector([cos(a1(1)); sin(a1(1)); -a2],Chh);
                
            case 2
                %% (2) Euclidean coordinate vector
                % homogeneous coordinates
                [h,Chh] =   sugr_Line_2D_Hes2hom(a1,a2);
                % minimal representation
                Line_2D   = sugr_minimal_vector(h,Chh);
                
            case 3
                %% (3) homogeneous coordinates
                % minimal representation
                Line_2D=sugr_minimal_vector(a1,a2);
                
        end
    case 3
        size_a3=size(a3,1);
        switch size_a3
            case 1
                %% (3) homogeneous coordinates, CovM == 0
                % minimal representation
                Line_2D = sugr_minimal_vector([a1;a2;a3],Chh);
                
            case 2
                %% (3) Euclidean coordinates
                % spherically normalized homogeneous vector
                [h,Chh] = sugr_Line_2D_Hes2hom([a1;a2],a3);
                % minimal representation
                Line_2D = sugr_minimal_vector(h,Chh);
                
        end
        
    case 4
        %% (4) homogeneous coordinates
        % minimal representation
        Line_2D = sugr_minimal_vector([a1;a2;a3],a4);
        
    case 5
        %% (5) centroid form
        % homogeneous coordinates
        [h,Chh] = sugr_Line_2D_cen2hom([a1;a2],a3,a4,a5);
        Line_2D = sugr_minimal_vector(h,Chh);
end
Line_2D.Jr = null(Line_2D.h');
Line_2D.type = 2;
