%% Minimal representation of Projection P 3D-->2D
%
% P = sugr_minimal_Projection_3D_2D(P,Cpp)
%
% Wolfgang Förstner 2/2013
% wfoerstn@uni-bonn.de
%
% See also sugr_Projection_3D_2D, sugr_minimal_vector

function Ps = sugr_minimal_Point_Projection_3D_2D(P,Cpp)

p = P(:);                                  % vectorize
p_minimal  = sugr_minimal_vector(p,Cpp);   % determine minimal vector
Ps.P       = reshape(p_minimal.h,3,4);     % build P-matrix
Ps.Crr     = p_minimal.Crr;                % store CovM_rr
