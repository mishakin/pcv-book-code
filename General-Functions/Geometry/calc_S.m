%% skew matrix of 3-vector
% see PCV (A.34)
% 
% Usage:
%   S = calc_S(x)
%
%   x - 3x1 vector
%   S - 3x3 skew matrix
%
% let x = [a b c], than
%     S = [  0 -c  b |
%         | -c  0 -a |
%         | -b  a  0 ]
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de 
%
% See also calc_S_reduced, calc_v_from_S

function S = calc_S(x)
S = [...
     0    -x(3)  x(2);...
     x(3)    0  -x(1);...
    -x(2)  x(1)    0 ...
    ];
