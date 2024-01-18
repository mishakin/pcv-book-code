%% K, R, Z, K33 from P (without covariance propagation)
% see Algorithm 18, PCV, p. 500, without covariance propagation
%
% Usage:
%   [K,R,Z,K33] = calc_KRZ_from_P(P,sc)
%  
%   P   - 4 x 3-projection matrix P=K*R*[I|-Z]
%   sc  - sign of the principle distance +1, -1, optional, default sc = -1;
% 
%   K   - calibration matrix 
%   R   - rotation matrix
%   Z   - projection centre
%   K33 - sign of K33 from the decomposition
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% See also calc_P_from_KRZ, calc_Q_from_P, calc_viewing_direction

function [K,R,Z,K33] = calc_KRZ_from_P(P,sc)

if nargin == 1
    sc = -1;
end

% partitioning
A = P(1:3,1:3);
a = P(1:3,4);
Ainv = inv(A);

% projection centre
Z = -Ainv*a;

% normalize P
sdetA = sign(det(A));
A = A*sdetA;

% calibration and rotation matrix
[Rinv,Kinv] = qr(inv(A));
R  = Rinv';
K0 = inv(Kinv);

% change chirality abd signs in K if necessary
Di = diag(sign(diag(K0)))*diag([sc,sc,1]);
R  = Di*R;
K0 = K0*Di;                                                                %#ok<MINV>

% nomalize K
K33 = K0(3,3); %*sign(sc);
K   = K0/K33;
