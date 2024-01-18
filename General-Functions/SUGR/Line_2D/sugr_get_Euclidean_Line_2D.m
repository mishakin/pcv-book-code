%% Determines Hessian paramters of 2D line
%
% [e,Cee] = sugr_get_Euclidean_Line_2D(l)
%
% * l = 2D line struct, minimal representation {l.h,l.Crr}
%
% * e = 2-vector with direction of normal and distance (Hessian)
% * Cee = CovM of e
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de
%
% See also sugr_Line_2D, sugr_Line_2D_hom2Hes, sugr_Line_2D_hom2cen,
% sugr_Line_2D_Hes2hom, sugr_Line_2D_Hes2cen, sugr_Line_2D_cen2hom,
% sugr_get_centroid_Line_2D

function [e,Cee] = sugr_get_Euclidean_Line_2D(l)

switch sugr_get_isfinite_Line_2D(l)
    case 0  % infinite line
        e   = zeros(2,1);
        Cee = zeros(2);
    case 1  % finite line
        n    = norm(l.h(1:2));                                     % norm
        e    = [atan2(l.h(2),l.h(1)),-l.h(3)/n]';                  % Euclidean parameters
        J    = 1 / n^3 * ...                                       % J: Jeh * Jhr
            [     -n * l.h(2) ,       n * l.h(1) ,  0;   ...
            l.h(3) * l.h(1) ,  l.h(3) * l.h(2) , -n^2] ...
            * l.Jr;
        Cee  = J * l.Crr * J';
end