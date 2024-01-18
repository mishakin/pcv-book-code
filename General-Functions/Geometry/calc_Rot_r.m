%% Cayley rotation matrix
% R = (I-S_r)^-1*(I+S_r)
% see PCV (8.62)
% 
% Usage:
%   R = calc_Rot_r(r)
% 
%   r - double 3x1
%   R - double 3x3, 3D rotation matrix (rotates by 2|r|)
%
% Wolfgang Förstner 10/2011
% wfoerstn@uni-bonn.de 
%
% See also calc_R_from_opq_kraus, calc_r_phi_from_R, calc_Rot_ab,
% calc_Rot_q, calc_Rot_rod, calc_Rot_rp, calc_Mq, calc_Mq_comm, calc_opq_from_R_Kraus

function R = calc_Rot_r(r)

R = inv(eye(3)-calc_S(r))*(eye(3)+calc_S(r));                               %#ok<MINV>
