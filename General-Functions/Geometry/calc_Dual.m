%% dualizing matrix
%
% Usage:
%   D = calc_Dual
%
%   D - 6x6 dualizing matrix
%       [0 0 0 1 0 0
%        0 0 0 0 1 0
%        0 0 0 0 0 1
%        1 0 0 0 0 0
%        0 1 0 0 0 0
%        0 0 1 0 0 0]
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% See also calc_Gamma, calc_Gammadual, calc_Gamma_reduced, 
% calc_Gammadual_reduced, calc_Pi calc_Pidual

function D = calc_Dual

D = [zeros(3,3),eye(3); 
     eye(3),zeros(3,3)
    ];

