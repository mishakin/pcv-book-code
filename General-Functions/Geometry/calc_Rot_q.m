%% Rotation matrix from quaternion 
% q = [q_0,qvector]
% see PCV (8.55)
%
% Usage:
%   R = calc_Rot_q(q)
% 
%   q - double 4x1, quaternion
%   R - double 3x3, 3D rotation matrix
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% See also calc_R_from_opq_kraus, calc_r_phi_from_R, calc_Rot_ab,
% calc_Rot_r, calc_Rot_rod, calc_Rot_rp, calc_Mq, calc_Mq_comm, calc_opq_from_R_Kraus

function R = calc_Rot_q(q)

q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

R = [
  q0^2+q1^2-q2^2-q3^2    2*(q1*q2-q3*q0)    2*(q1*q3+q2*q0); ... 
    2*(q1*q2+q3*q0)    q0^2-q1^2+q2^2-q3^2  2*(q2*q3-q1*q0); ...
    2*(q1*q3-q2*q0)      2*(q2*q3+q1*q0)   q0^2-q1^2-q2^2+q3^2 ...
    ] ...
    /(q0^2+q1^2+q2^2+q3^2);
