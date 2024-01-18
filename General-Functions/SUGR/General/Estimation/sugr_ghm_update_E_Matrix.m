%% GHM update params relative orientation/E-matrix
%
% estx = sugr_ghm_update_E_Matrix(estx,estx_r);
%
% estx = 3 x 4 matrix [b,R]
% estxr = 1 x 5 vector [db_r;dr]
%
% Wolfgang Förstner 09/2011
% wfoerstn@uni-bonn.de 

function estx = sugr_ghm_update_E_Matrix(estx,estx_r)

% initial values
ba = estx(:,1);
Ra = estx(:,2:4);

% updates
best = ba + null(ba') * estx_r(1:2);
best = best / norm(best);
Rest = expm(calc_S(estx_r(3:5)')) * Ra;

% updated parameters
estx = [best, Rest];
