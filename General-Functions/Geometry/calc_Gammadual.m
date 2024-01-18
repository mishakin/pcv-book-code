%% dual Plückermatrix from Plücker vector
% see PCV (5.134)
%
% Usage
%   G = calc_Gamma (L)
%
%   L - 6x1 Plücker Line
%   G - 4x1 dual Gamma Matrix / dual Plückermatrix
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% See also calc_Gamma, calc_Gamma_reduced, calc_Gammadual_reduced, 
% calc_Pi, calc_Pidual, calc_Dual


function G = calc_Gammadual (L)

G= ...
[...
       0     L(3) -L(2) -L(4); ...
      -L(3)  0     L(1) -L(5); ...
       L(2) -L(1)  0    -L(6); ...
       L(4)  L(5)  L(6)  0 ...
];

