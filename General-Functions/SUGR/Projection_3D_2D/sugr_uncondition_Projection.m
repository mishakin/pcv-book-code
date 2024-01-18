%% Uncondition Projection
%
% P = sugr_uncondition_Projection(Pc, Ml, Mr) 
%
% Pc = struct, projection conditioned
%      Pc.P   3 x 4-matrix
%      Pc.Crr reduced CovM
% Ml = 3 x 3 condition matrix
% Mr = 4 x 4 condition matrix
% 
% P  = struct, unconditioned matrix P = inv(Ml) * Pc * Mr, 
%
% Wolfgang Förstner 2/2013
% wfoerstn@uni-bonn.de
%
% See also sugr_condition_Points

function P = sugr_uncondition_Projection(Pc, Ml, Mr)

% Matrix
Pmc = Pc.P;

% unconditioned matrix
Pm =  inv(Ml) * Pmc * Mr;                                                  %#ok<MINV>
% Pm_unconditioned = Pm;

% uncondition CovM
pc  = Pmc(:);
p   = Pm(:);
% dpu    =                    (Mr' kron inv(Ml)) * dpc
% dp_r   =          J'(dp') * (Mr' kron inv(Ml)) * J(pc') dpc_r
% dp_h   = J(dp') * J'(dp') * (Mr' kron inv(Ml)) * J(pc') dpc_r
J = null(p') * null(p')' * kron(Mr',inv(Ml)) * null(pc');
Chh = J * Pc.Crr * J';

P = sugr_Projection_3D_2D(Pm,Chh);


