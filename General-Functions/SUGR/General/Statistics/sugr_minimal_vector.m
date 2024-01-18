%% generate minimal representation for homogeneous vector
%
% xs = sugr_minimal_vector(x,Cxx)
%
% * x    uncertain homogeneous vector
% * Cxx  covariance matrix
% * xs =   
% * .h     spherically normalized
% * .Crr   reduced CovM
% * .Jr    Jacobian xr -> xs 
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de 
%
% See also sugr_minimal_3D_Line

function xs = sugr_minimal_vector(x,Cxx)

n      = norm(x);             % norm
xs.h   = x/n;                 % spherical normalization
Jr     = null(xs.h');         % Jacobian xr -> xs 
Crr    = Jr' * Cxx * Jr/n^2;  % CovM of reduced vector
xs.Crr = Crr;                 % Crr
xs.Jr  = Jr;  

