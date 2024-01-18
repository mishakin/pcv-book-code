%% Plot 2D point with standard ellipse
%
% sugr_plot_Point_2D(x,center_type,bound_type,lw,factor)
%
% * x      = geometric element
% * center_type   = type of plotting the element
% * bound_type    = type of plotting the confidence region
% * lw      = linewidth
% * factor  = magnification factor
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de
%
% sw 9/2016
%
% See also sugr_Point_2D

function sugr_plot_Point_2D(x,center_type,bound_type,lw,factor)


if nargin < 4
    lw=2;
end
if nargin < 5
    factor=1;
end

if sugr_get_isfinite_Point_2D(x) % only finite elements
    Chh = sugr_get_CovM_homogeneous_Vector(x);
    C   = adjunctMatrix(Chh-x.h*x.h');
    sugr_plot_conic_explicit(C,center_type,bound_type,lw,factor);
    
end


