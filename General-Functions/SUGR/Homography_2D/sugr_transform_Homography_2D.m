%% transform uncertain homography with crisp transformations
% 
% H = sugr_transform_Homography_2D(H0,Tl,Tr);
%
% H = Tl * H0 * Tr
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 3/2011

function H = sugr_transform_Homography_2D(H0,Tl,Tr)

% homography
Ht = Tl * H0.H * Tr;

% covariance matrix
C0hh = sugr_get_CovM_homogeneous_Homography_2D(H0);
J = kron(Tr', Tl);
Chh = J * C0hh *J';

% generate struct
H = sugr_Homography_2D(Ht,Chh);
