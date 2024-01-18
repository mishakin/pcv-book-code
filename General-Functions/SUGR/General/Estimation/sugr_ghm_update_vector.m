%% GHM update spherically normalized vector
%
% xu = sugr_ghm_update_vector(x,p)
%
% * x = N-vector, approximate values
% * p = (N-1)-vector, reduced corrections
%
% xu = updated vector
%
% Wolfgang Förstner 03/2011
% wfoerstn@uni-bonn.de 

function xu = sugr_ghm_update_vector(x,p)

% Update in tangent space
xu = x + null(x') * p;

% normalize
xu = xu / norm(xu);
