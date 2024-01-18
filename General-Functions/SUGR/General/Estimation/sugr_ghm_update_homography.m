%% GHM update spectrally normalized homography
%
% Hu = sugr_ghm_update_homography(Ha,k)
%
% * Ha = N x N matrix, approximate values
% * k  = N^2-1-vector, reduced corrections 
%
% Hu = updated matrix
%
% Wolfgang Förstner 03/2011
% wfoerstn@uni-bonn.de 

function Hu = sugr_ghm_update_homography(Ha,k)

% dimension
N = size(Ha,1);

% K = traceless matrix
K      = reshape([k;0],N,N);
K(N,N) = -trace(K);

% update
Hu = expm(K)*Ha;

