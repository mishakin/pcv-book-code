%% calculates P from K, R and Z
% see PCV (12.34), and additional sign of K33
%
% Usage:
%   P = calc_P_from_KRZ (K,R,Z,[K33])
%  
%   K   - calibration matrix 
%   R   - rotation matrix
%   Z   - projection centre
%   K33 - sign of K33, optional, default K33 = 1;
%
%   P   - 4 x 3-projection matrix P = K*R*[I|-Z] 
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% See also calc_KRZ_from_P, calc_Q_from_P, calc_viewing_direction

function P = calc_P_from_KRZ (K,R,Z,K33)

if nargin ==3
    K33=1;
end

P = K*R*[eye(3) -Z]*K33;
