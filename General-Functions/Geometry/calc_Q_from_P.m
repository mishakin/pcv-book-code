%% calculates projection matrix Q for lines from P
% see PCV (12.71)
%
% Usage:
%   Q = calc_Q_from_P (P)
%  
%   P - 4 x 3-projection matrix P = K*R*[I|-Z] 
%   Q - 3x6 projection matrix for lines
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% See also calc_P_from_KRZ, calc_KRZ_from_P, calc_viewing_direction

function Q = calc_Q_from_P (P)

A=P(1,:)';
B=P(2,:)';
C=P(3,:)';

Q=[(calc_Pi(B) * C)'; ...
   (calc_Pi(C) * A)'; ...
   (calc_Pi(A) * B)'];
