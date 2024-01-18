%% cofactor of a square matrix
% 
% B = cofactorMatrix(A)
%
% inv(A) = cof(A)'/det(A)   holds.
%   d det(A)/ dA = cof(A)   holds.
%
% author:  J. Meidow, FGAN-FOM
%
% $Log: cof.m,v $
% Revision 1.1  2007/09/26 06:33:52  meidow
% *** empty log message ***
%
% adapted Wolfgang Foerstner 1/2011
% wfoerstn@uni-bonn.de 

function B = cofactorMatrix(A)

B = adjunctMatrix(A)';

