%% Uncertain centered line parameters from uncertain Hessian line parameters hes
%
% [x0,p,sp,sq]= sugr_Line_2D_Hes2cen(e,Cee)
%
% * e    = line parameters [phi, d]'
% * Cee  = 2 x 2 covariance matrix
%
% * x0 = centroid, [x-coord. y-coord.]'
% * p  = direction of normal
% * sp = standard deviation of normal direction
% * sq = standard deviation across line
%
% Jochen Meidow, Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de 
%
% See also sugr_Line_2D, sugr_Line_2D_hom2Hes, sugr_Line_2D_hom2cen,
% sugr_Line_2D_Hes2hom, sugr_Line_2D_cen2hom, 
% sugr_get_Euclidean_Line_2D, sugr_get_centroid_Line_2D

function [x0,p,sp,sq]= sugr_Line_2D_Hes2cen(e,Cee)

p   = e(1);                         % phi
d   = e(2);                         % d

R = [cos(p),-sin(p) ;...
     sin(p), cos(p)];               % rotation 
if Cee(1,1) ~= 0
    u0 = Cee(1,2)/Cee(1,1);         % offset
else 
    u0 = 0;
end

xy    = R*[d;u0];                   % x0, y0

x0 = [xy(1),xy(2)]';                % vector x0
sp = sqrt(Cee(1,1));                % sigma phi
sq = sqrt(-u0^2*Cee(1,1)+Cee(2,2)); % sigma q

