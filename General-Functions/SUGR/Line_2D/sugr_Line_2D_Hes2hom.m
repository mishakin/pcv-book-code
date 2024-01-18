%% Uncertain homogenous line parameters from uncertain Hessian (Euclidean) line parameters
%
% [h,Chh] = sugr_Line_2D_Hes2hom(e,Cee) 
%
% * e    = line parameters [phi, d]'
% * Cee  = 2 x 2 covariance matrix
% * h    =  homogeneous line vector, spherically normalized
% * Chh  = 3 x 3 covariance matrix
%
% Wolfgang Förstner  1/2011
% wfoerstn@uni-bonn.de 
%
% See also sugr_Line_2D, sugr_Line_2D_hom2Hes, sugr_Line_2D_hom2cen,
% sugr_Line_2D_Hes2cen, sugr_Line_2D_cen2hom, 
% sugr_get_Euclidean_Line_2D, sugr_get_centroid_Line_2D

function [h,Chh] = sugr_Line_2D_Hes2hom(e,Cee)


h    = [cos(e(1)); ...  
        sin(e(1)); ...
           -e(2)];                  % homogeneous vector
J    = [-sin(e(1)) ,    0;...       % Jacobian 
         cos(e(1)) ,    0;...
                   0 ,   -1 ];
Chh = J * Cee * J';                 % covariance matrix  

