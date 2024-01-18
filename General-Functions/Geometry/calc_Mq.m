%% Multiplication matrix for quaternions
% q p = M(q)*p, see PCV (8.45)
%
% Usage:
%	M = calc_Mq(q)
%
%   q - double 4x1, quaternion
%
%   M - double 4x4, Multiplication matrix for quaternions
%
% Wolfgang F�rstner
% wfoerstn@uni-bonn.de 
%
% See also calc_Mq_comm, calc_R_from_opq_kraus, calc_r_phi_from_R, 
% calc_Rot_ab, calc_Rot_q, calc_Rot_r, calc_Rot_rod, calc_Rot_rp, 
% calc_opq_from_R_Kraus

function M = calc_Mq(q)

a = q(1);
b = q(2);
c = q(3);
d = q(4);

M=[ a -b -c -d; ...
    b  a -d  c; ...
    c  d  a -b; ...
    d -c  b  a
    ];