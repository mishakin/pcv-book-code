%% Plot uncertain 2D line with standard hyperbola
%
% sugr_plot_Line_2D(l,center_type,bound_type,lw,factor,type)
%
% * l      = geometric element
% * center_type   = type of plotting the element
% * bound_type    = type of plotting the confidence region
% * lw      = linewidth
% * factor  = magnification factor
% * se      = [s,e]'  start/end for line segment
%             [0,0]'  then line segment (+-sqrt(2))
%             [1,1]'  then line until xlim,ylim (default)
%
% Wolfgang Förstner  1/2011
% wfoerstn@uni-bonn.de
%

function sugr_plot_Line_2D(l,center_type,bound_type,lw,factor,se)

if nargin < 6
    se=[1,1];
end

if sugr_get_isfinite_Line_2D(l)     % only finite elements
    
    if nargin < 4
        lw=2;
    end
    
    if nargin < 5
        factor=1;
    end
    
    Chh = sugr_get_CovM_homogeneous_Vector(l);
    C   = Chh-l.h*l.h';
    
    [~,phi,~,~] = sugr_get_centroid_Line_2D(l);
    
    if phi < pi/4 && phi > -3*pi/4 && nargin == 6
        se = -se;
    end
    
    sugr_plot_conic_explicit(C,center_type,bound_type,lw,factor,se);
    
end
