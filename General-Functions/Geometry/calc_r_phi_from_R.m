%% rotation axis and rotation angle from rotation matrix R
% see PCV Sect. 8.1.4.2, p. 331 ff
%
% Usage:
%   [r,phi]= calc_r_phi_from_R(R)
%
%   R - double 3x3, 3D rotation matrix
%
%   r - double 3x1, normalized rotation axis
%   phi - double, rotation angle, radiant
%
% Wolfgang Förstner 10/2011
% wfoerstn@uni-bonn.de 
%
% See also calc_R_from_opq_kraus, calc_opq_from_R_Kraus, calc_Rot_ab, calc_Rot_q, 
% calc_Rot_r, calc_Rot_rod, calc_Rot_rp, calc_Mq, calc_Mq_comm

function [r,p]= calc_r_phi_from_R(R)

tr = trace(R);

% vector a which will lead to the rotation axis below, PCV Eq.(8.33)
a = -[R(2,3)-R(3,2),R(3,1)-R(1,3),R(1,2)-R(2,1)]';
an = norm(a);

% rotation angle, PCV Eq.(8.34)
phi = atan2(an,tr-1);

% three cases of rotation axis
if abs(phi) < eps 
   % 0°
   r = [1,0,0]';
   p = 0;
elseif abs(phi-pi) < eps
   % 180°
   D = (R+eye(3));       
   [~,i] = max([D(1,1),D(2,2),D(3,3)]);
   ru = D(:,i);
   r = ru/norm(ru);
   p = pi;
else
   % general case 
   r = a/an;
   p = phi;
end

