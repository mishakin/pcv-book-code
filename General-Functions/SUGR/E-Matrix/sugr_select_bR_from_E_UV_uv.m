%% Select bR from E
% select one of the four cases;
%
% similar to Alg. 20, but from (U,V) instead from E
%
% [b,R,zwe] = sugr_select_bR_from_E(U0,V0,xl,xr);
%
% U0, V0    = svd(E), with positive determinants enforced
% xl,xr     = Nx3 matrices for corresponding points
%
% Wolfgang Förstner 9/2011
% wfoerstn@uni-bonn.de
%
% See also sugr_E_Matrix

function [be,Re] = sugr_select_bR_from_E_UV_uv(U,V,u,v)

% auxiliary matrices
W  = [0,1,0;-1,0,0;0,0,1];
u3 = U(:,3);

% number of corresponding directions
N = size(u,1);

% initialize 
solution = 0;

be       = [1,0,0]';
Re       = eye(3);

% determine alternatives (PCV p. 583, Alg. 20) 
for t=0:3
    switch t
%         case 0
%             R  = V*W*U';     % R'=U*W*V' (HZ) -> R = V*W'*U'
%             b  = u3;
%         case 1
%             R   = V*W*U';
%             b  = -u3;
%         case 2
%             R   = V*W'*U';
%             b  = u3;
%         case 3
%             R   = V*W'*U';
%             b  = -u3;
    case 0
        R  = V*W'*U';     % R'=U*W*V' (HZ) -> R = V*W'*U'
        b  = u3;
    case 1
        R   = V*W'*U';   
        b  = -u3;
    case 2
        R   = V*W*U';
        b  = u3;
    case 3
        R   = V*W*U';
        b  = -u3;
    end   
    S = calc_S(b);
%% check for positiveness (not yet included: check for direction || basis )
    % m-vectors: 3 x N Matrix, PCV (13.165)
    Nm = -S*S*u(:,1:3)';   
    % m x (u x R'v) - vectors, 3 x M-matrix, part of denom., PCV (13.164) 
    Mm = (calc_S( Nm(1:3,:)) * cross(u(:,1:3)', R'*v(:,1:3)'))';
    % 
    ssn = sign((ones(N,1) * b' .* Mm) * [1 1 1]');
    srn = ssn.* sign(((R'*v(:,1:3)')' .* Nm(1:3,:)')*[1 1 1]'); 
    if mean(srn) > 0 && mean(ssn) > 0 
        solution = 1;
        be       = calc_v_from_S(S);
        Re       = R;
    end
end
if solution == 0
    be = zeros(3,1);
end
