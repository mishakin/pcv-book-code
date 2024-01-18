%% Test whether 2D line is finite
%
% f = sugr_get_isfinite_Line_2D(l)
%
% * l    = 2D line struct {l.h,l.Crr}
% * f    = boolean, true if l is finite, false otherwise
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de 
%
% See also sugr_Line_2D

function f = sugr_get_isfinite_Line_2D(l)

global Threshold_Euclidean_Normalization

f   =   norm(l.h(1:2)) > ...
                Threshold_Euclidean_Normalization * abs(l.h(3));
            
