%% calculate Euler angles from Rotation matrix 
% following PCV Sect. 8.1.3.3 (8.23) ff
%
% Usage
%   [omega,phi,kappa] = calc_R_from_opq_kraus(R)
% 
%   R - double 3x3, 3D rotation matrix
% 
%   omega,phi,kappa - Euler angles, radiant
%
% R = R_3(kappa)*R_2(phi)*R_1(omega)
%   = [...
%     cos(phi)*cos(kappa), ...
%     -cos(phi)*sin(kappa), ...
%     sin(phi); ... 
%     sin(omega)*sin(phi)*cos(kappa)+cos(omega)*sin(kappa), ...
%     -sin(omega)*sin(phi)*sin(kappa)+cos(omega)*cos(kappa), ...
%     -sin(omega)*cos(phi); ...
%     -cos(omega)*sin(phi)*cos(kappa)+sin(omega)*sin(kappa), ...
%     cos(omega)*sin(phi)*sin(kappa)+sin(omega)*cos(kappa), ...
%     cos(omega)*cos(phi)];
%
% o      = omega <= pi
% p      = phi
% q      = kappa
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% See also calc_R_from_opq_kraus, calc_r_phi_from_R, calc_Rot_ab, calc_Rot_q, 
% calc_Rot_r, calc_Rot_rod, calc_Rot_rp, calc_Mq, calc_Mq_comm

function [o,p,q] = calc_opq_from_R_Kraus(R)

q = mod(atan2(-R(1,2),R(1,1))+2*pi,2*pi);
o = mod(atan2(-R(2,3),R(3,3))+2*pi,2*pi);
p = mod(atan2(R(1,3),-R(1,2)/sin(q))+2*pi,2*pi);

if  o > pi
    q = q-pi;
    o = o-pi;
    p = mod(pi-p+2*pi,2*pi);
end