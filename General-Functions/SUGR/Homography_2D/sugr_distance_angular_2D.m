%% Angular distance of two 2D entities
%
% distance = sugr_distance_angular_2D(a,b)
%
% *       'Point_2D'+'Point_2D'
% *       'Point_2D'+'Line_2D'
% *       'Line_2D'+'Point_2D'
%
% TODO: sugr_distance_Euclidean_2D 'Line_2D'+'Line_2D' if parallel
%
% distance.a  = angular distance
%
% distance.sa = standard deviation of angular distance
%
% distance.dof = degrees of freedom = 1
%
% if distance = 0 -> dof =-1
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 2/2011

function distance = sugr_distance_angular_2D(a, b)

global Threshold_Euclidean_Normalization

pair = strcat(num2str(a.type), num2str(b.type));

distance.dof = 1;

if norm(a.h - b.h) > Threshold_Euclidean_Normalization
    switch pair
        case '11'
            Cahh = sugr_get_CovM_homogeneous_Vector(a);
            Cbhh = sugr_get_CovM_homogeneous_Vector(b);
        case '12'
            Cahh = sugr_get_CovM_homogeneous_Vector(a);
            Cbhh = sugr_get_CovM_homogeneous_Vector(b);
        case '21'
            Cahh = sugr_get_CovM_homogeneous_Vector(a);
            Cbhh = sugr_get_CovM_homogeneous_Vector(b);
        case '22'
            Cahh = sugr_get_CovM_homogeneous_Vector(a);
            Cbhh = sugr_get_CovM_homogeneous_Vector(b);
    end
    ncab = norm(cross(a.h, b.h));         % norm crossproduct
    ipab = a.h'*b.h;                      % inner product
 
    angle = atan2(ncab, ipab);            % angle
    distance.a = angle;
 
    ja = (ipab * a.h'-b.h') / ncab; % Jacobian angle --> a.h
    jb = (ipab * b.h'-a.h') / ncab; % Jacobian angle --> b.h
 
    distance.sa = sqrt(ja * Cahh * ja' + jb * Cbhh * jb');
else
    distance.a = 0;
    distance.sa = 0;
    distance.dof = - 1;
end

