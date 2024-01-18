%% Get Euclidean coordinates with CovM
%
% [e,Cee] = sugr_get_Euclidean_Point_2D(x)
% 
% if finite point:
%    e   = Euclidean coordinates
%    Cee = CovM of e
% else zeros
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
% wf 1/2011
%
% See also sugr_Point_2D, sugr_get_Euclidean_Line_2D,
% sugr_get_Euclidean_Point_3D, sugr_get_isfinite_Point_2D

function [e,Cee] = sugr_get_Euclidean_Point_2D(x)

switch sugr_get_isfinite_Point_2D(x)
    case 0
        e   = zeros(2,1);
        Cee = zeros(2,2);
    case 1
        e   = x.h(1:2)/x.h(3);
        J   =  1/x.h(3) * [eye(2) -x.h(1:2)/x.h(3)] * null(x.h');
        Cee = J * squeeze(x.Crr)  * J';
end
