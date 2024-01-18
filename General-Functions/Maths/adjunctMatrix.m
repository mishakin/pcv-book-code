%% determines adjunct matrix for a square matrix
% 
% B = adjunctMatrix(A)
%
% inv(A) = adj(A)/det(A)   holds.
%   d det(A)/ dA = adj(A)'   holds.
%
% author:  J. Meidow, FGAN-FOM
%
% $Log: adj.m,v $
% Revision 1.1  2007/09/26 06:33:52  meidow
% *** empty log message ***
%
% adapted Wolfgang Foerstner 1/2011
% wfoerstn@uni-bonn.de 

function B = adjunctMatrix(A)

B = zeros(size(A));
for i=1:size(A,1)
    idx_i = setdiff( 1:size(A,1), i);
    for j=1:size(A,2)        
        idx_j = setdiff( 1:size(A,2),j);
        B(i,j) = (-1)^(i+j)*det(A(idx_j,idx_i));
    end
end
