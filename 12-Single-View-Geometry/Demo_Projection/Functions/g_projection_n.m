%% Check the Jacobian of the projection constraint
%
% g  = N(x') - N(P X)
%
% px = [vecP; X; x']
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 04/18
% wenzel@igg.uni-bonn.de
%
function g = g_projection_n(px)

yn = px(17:19); 
yn = yn/norm(yn);

ypn = reshape( px(1:12),3,4 ) * px(13:16);
ypn = ypn/norm(ypn);

g = yn - ypn;

