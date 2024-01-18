%% reduced dual Plückermatrix from Plücker vector
% see PCV (7.132) ff
%
% Usage
%   [Gr , Rm] = calc_Gammadual_reduced(L)
%
%   L - 6x1 Plücker Line
%   Gr - 2x4 reduced dual Gamma matrix / reduced dual Plückermatrix
%   Rm - 2x4 reduction matrix
%
% Wolfgang Förstner 02/2010
% wfoerstn@uni-bonn.de 
%
% See also calc_Gamma, calc_Gammadual, calc_Gamma_reduced, 
% calc_Pi, calc_Pidual, calc_Dual

function [Gr , Rm] = calc_Gammadual_reduced(L)

[~,i] = max(abs(L));
G = calc_Gammadual(L);

switch i 
    case 1
        Rm = [0 1 0 0;0 0 1 0];
        Gr = G([2,3]',:);
    case 2
        Rm = [0 0 1 0;1 0 0 0];
        Gr = G([3,1]',:);
    case 3
        Rm = [1 0 0 0;0 1 0 0];
        Gr = G([1,2]',:);
    case 4
        Rm = [1 0 0 0;0 0 0 1];
        Gr = G([1,4]',:);
    case 5
        Rm = [0 1 0 0;0 0 0 1];
        Gr = G([2,4]',:);
    case 6
        Rm = [0 0 1 0;0 0 0 1];
        Gr = G([3,4]',:);
end
