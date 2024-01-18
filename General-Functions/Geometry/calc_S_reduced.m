%% Skew matrix with independent selected rows
% see PCV (7.120) ff
%
% Usage:
%   [S,R] = calc_S_reduced(x)
%
%   Ss - 2x3 skew matrix, reduced to independent rows
%   Rm - 2x3 reduction matrix: Ss = Rm * S 
%
% Wolfgang Förstner 02/2010
% wfoerstn@uni-bonn.de 
%
% See also calc_S, calc_v_from_S

function [S , R] = calc_S_reduced(x)

% select index with maximum |x|
[~,i] = max(abs(x));

% generate rselection matrix
switch i
    case 1
        R = [0 1 0;0 0 1];
    case 2
        R = [1 0 0;0 0 1];
    case 3
        R = [1 0 0;0 1 0];
end
% generate selected matrix
S = R *calc_S(x);

