%% 3-vector from skew matrix S
% see PCV (A.34)
%
% Usage:
%   x = calc_v_from_S(S)
%
%   S - 3x3 skew matrix
%   x - 3x1 vector
%
% Wolfgang Förstner wf 1/2011
% wfoerstn@uni-bonn.de 
%
% See also calc_S, calc_S_reduced

function x = calc_v_from_S(S)
    x=[S(3,2); S(1,3); S(2,1)];
