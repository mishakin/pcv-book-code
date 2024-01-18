%% Uncertain homogenous line parameters from uncertain centered line parameters
%
% [h,Chh] = sugr_Line_2D_cen2hom(x0,p,sp,sq)
% 
% * x = centroid, [x-coord. y-coord.]'
% * p  = direction of normal
% * sp = standard deviation of normal direction
% * sq = standard deviation across line
%
% * h = homogeneous vector
% * Chh = CovM
%
% Wolfgang Förstner  1/2011
% wfoerstn@uni-bonn.de 
%
% See also sugr_Line_2D, sugr_Line_2D_hom2Hes, sugr_Line_2D_hom2cen,
% sugr_Line_2D_Hes2hom, sugr_Line_2D_Hes2cen, 
% sugr_get_Euclidean_Line_2D, sugr_get_centroid_Line_2D

function [h,Chh] = sugr_Line_2D_cen2hom(x0,p,sp,sq)

ly  = [1,0,0]';                          % negative y-axis
Cyy = diag([0,sp^2,sq^2]);               % standardized Cov.-Matrix

TR  = inv([cos(p),-sin(p),x0(1);...
           sin(p), cos(p),x0(2);...
           0     , 0     , 1   ]');      % Transformation-matrix
h   = TR * ly;                           %#ok<*MINV> % hom. line
Chh = TR * Cyy * TR';                    % error-prop.

