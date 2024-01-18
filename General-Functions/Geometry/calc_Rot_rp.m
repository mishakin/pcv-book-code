%% Rotation matrix from axis and angle
% see PCV (8.29)
% 
% Usage:
%   R = calc_Rot_rp(r,p)
% 
%   r - double 3x1, rotation axis, normalized, such that |r| = 1
%   p - rotation angle [radiant]
%
%   R - double 3x3, 3D rotation matrix
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% See also calc_R_from_opq_kraus, calc_r_phi_from_R, calc_Rot_ab,
% calc_Rot_q, calc_Rot_r, calc_Rot_rod, calc_Mq, calc_Mq_comm, calc_opq_from_R_Kraus

function R = calc_Rot_rp(r,p)

rn = r/norm(r);
R = cos(p)*eye(3) + sin(p)*calc_S(rn) + (1-cos(p))*rn*rn';

