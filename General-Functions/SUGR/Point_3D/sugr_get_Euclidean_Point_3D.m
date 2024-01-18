%% Get Euclidean coordinates of 3D point with CovM
%
% [e,Cee] = sugr_get_Euclidean_Point_3D(x)
%
% if finite point:
%      e   = 3x1 Euclidean coordinates
%      Cee = 3x3 CovM of e
%
% else zeros
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de
%
% See also sugr_Point_3D

function [e,Cee] = sugr_get_Euclidean_Point_3D(x)

switch sugr_get_isfinite_Point_3D(x)
    case 0  % point at infinity
        e   = zeros(3,1);
        Cee = zeros(3);
    case 1  % finite point
        e   = x.h(1:3)/x.h(4);
        J   =  1/x.h(4) * [eye(3), -x.h(1:3)/x.h(4)] * null(x.h');
        Cee = J * squeeze(x.Crr)  * J';
end
