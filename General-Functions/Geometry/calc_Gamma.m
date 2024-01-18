%% Pl�ckermatrix from Pl�cker vector
% see PCV (5.72)
%
% Usage
%   G = calc_Gamma (L)
%
%   L - 6x1 Pl�cker Line
%   G - 4x1 Gamma Matrix / Pl�ckermatrix
%
% Wolfgang F�rstner
% wfoerstn@uni-bonn.de 
%
% See also calc_Gammadual, calc_Gamma_reduced, calc_Gammadual_reduced, 
% calc_P, calc_Pidual calc_Dual

function G = calc_Gamma (L)

G= ...
[...
       0     L(6) -L(5) -L(1); ...
      -L(6)  0     L(4) -L(2); ...
       L(5) -L(4)  0    -L(3); ...
       L(1)  L(2)  L(3)  0 ...
];

