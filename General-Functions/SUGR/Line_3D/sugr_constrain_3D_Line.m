%% imposes Plücker constraint for 3D lines
% 
% L = sugr_constrain_3D_Line(L) 
% 
% Wolfgang Förstner 7/2012
% wfoerstn@uni-bonn.de 
%
% See also sugr_Line_3D, sugr_minimal_3D_Line, sugr_construct_join_Line_3D

function L = sugr_constrain_3D_Line(L)

% take parts
Lh  = L(1:3);
L0  = L(4:6);

% possibly distorted directions
Dht = Lh/norm(Lh);
D0t = L0/norm(L0);

% auxiliary variables
d   = norm(Dht-D0t);
r   = sqrt(1-d^2/4);

% corrected disrections
Dh  = (1/2 + r/d)*Dht + (1/2-r/d)*D0t;
D0  = (1/2 - r/d)*Dht + (1/2+r/d)*D0t;

% corrected line parameters
L   = [norm(Lh)*Dh; norm(L0)*D0];


