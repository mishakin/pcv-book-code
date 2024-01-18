%% Test whether 3D point is finite
%
% f = sugr_get_isfinite_Point_3D(x)
%
% * x    = 3D point struct {x.h,x.Crr}
% * f    = boolean, true if x is finite, false otherwise
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de
%
% See also sugr_Point_3D

function f = sugr_get_isfinite_Point_3D(x)

global Threshold_Euclidean_Normalization

f   =  abs(x.h(4)) > ...
                Threshold_Euclidean_Normalization * norm(x.h(1:3));
            
