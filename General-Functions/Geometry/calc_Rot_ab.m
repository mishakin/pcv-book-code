%% minimum Rotation matrix from vector a to vector b
% see PCV (8.76)
%
% Usage:
%   R = calc_Rot_ab(a,b)
% 
%   a,b - double DX1, pair of noncoplanar vectors
%   R - double DxD, rotation matrix
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de 
%
% See also calc_R_from_opq_kraus, calc_opq_from_R_Kraus, calc_Rot_q, 
% calc_Rot_r, calc_Rot_rod, calc_Rot_rp, calc_Mq, calc_Mq_comm, calc_r_phi_from_R

function R = calc_Rot_ab(a,b)

as = a/norm(a);
bs = b/norm(b);

R = eye(length(a))+2*bs*as'-(as+bs)*(as+bs)'/(1+as'*bs); 