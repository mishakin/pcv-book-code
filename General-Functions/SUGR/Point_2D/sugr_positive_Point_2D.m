%% Change sign of 2D point if it is negative
% 
% pp = sugr_positive_Point2D(p)
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de
%
% See also sugr_Point_2D

function pp = sugr_positive_Point_2D(p)

pp = p;
if  sugr_get_isfinite_Point_2D(p)
    pp.h = p.h*sign(p.h(3));
end

    