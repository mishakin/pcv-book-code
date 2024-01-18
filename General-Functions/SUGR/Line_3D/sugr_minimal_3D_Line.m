%% generate minimal representation for 3D homogeneous line
%
% Ls = sugr_minimal_3D_Line(L,CLL)
%
% * L 6x1 uncertain homogeneous 3D line
% * CLL covariance matrix
%
% * Ls =   
% * .h spherically normalized
% * .Crr, reduced CovM
% * .Jr    Jacobian Lr -> Ls 
%
% Wolfgang Förstner 7/2012
% wfoerstn@uni-bonn.de 
%
% See also sugr_Line_3D, sugr_construct_join_Line_3D, sugr_constrain_3D_Line

function Ls = sugr_minimal_3D_Line(L,CLL)

% impose Plücker constraint, no change of length
L = sugr_constrain_3D_Line(L);

% generate 3D line structure, normalized spherically
Ls = sugr_minimal_vector(L,CLL);

% add Plücker constraint to CovM
h = calc_Dual*Ls.h; 
P = (eye(6) - h*h')/norm(L);    % scale to length 1
CLL = P * CLL * P;

% reduce covariance matrix
J6 =  null([Ls.h h]');
Ls.Crr = J6' * CLL * J6;
Ls.Jr  = J6;

