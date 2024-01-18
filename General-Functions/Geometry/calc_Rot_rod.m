%% Rotation matrix according to Rodriguez
% see PCV (8.60)
% 
% Usage:
%   R = calc_Rot_rod(q)
% 
%   q - double 3x1, Rodriguez representation
%   R - double 3x3, 3D rotation matrix
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% See also calc_R_from_opq_kraus, calc_r_phi_from_R, calc_Rot_ab,
% calc_Rot_q, calc_Rot_r, calc_Rot_rp, calc_Mq, calc_Mq_comm, calc_opq_from_R_Kraus

function R = calc_Rot_rod(q)

q1 = q(1)/2;
q2 = q(2)/2;
q3 = q(3)/2;

R = [
    1+q1^2-q2^2-q3^2  2*(q1*q2-q3)  2*(q1*q3+q2); ... 
    2*(q1*q2+q3)  1-q1^2+q2^2-q3^2  2*(q2*q3-q1); ...
    2*(q1*q3-q2)  2*(q2*q3+q1) 1-q1^2-q2^2+q3^2 ...
    ] ...
    /(1+q1^2+q2^2+q3^2);


