%% Check the Jacobian of the projection
% x' = N(PX)
%
% px = [vecP;X]
% y  = N(x')
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 04/18
% wenzel@igg.uni-bonn.de
%
function y = projection_n(px)

y = reshape( px(1:12),3,4 ) * px(13:16);
y = y/norm(y);

