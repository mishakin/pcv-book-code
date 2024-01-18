%% Pi-matrix from 3D-point 4-vector
% see PCV (7.37)
%
% Usage
%   P = calc_Pi (X)
%
%   X - 4x1 homogeneous 3D point
%   P - 6x4 Pi Matrix
%
% Wolfgang F�rstner
% wfoerstn@uni-bonn.de 
%
% See also calc_Pidual, calc_Gamma, calc_Gammadual, calc_Gamma_reduced, 
% calc_Gammadual_reduced, calc_Dual

function P = calc_Pi(X)
P= ...
[...
      X(4)  0     0    -X(1); ...
      0     X(4)  0    -X(2); ...
      0     0     X(4) -X(3); ...
      0    -X(3)  X(2)  0; ...
      X(3)  0    -X(1)  0; ...
     -X(2)  X(1)  0     0  ...
];
