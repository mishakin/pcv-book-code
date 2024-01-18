%% b and R from E using set of pairs of directions
%
% Alg. 20 
%
% [b, R, code] = sugr_bR_from_E_uv(E,u,v)
%
% E     = 3x3 E-matrix
% u,v   = Nx3 matrices of directions
%
% b     = 3x1 normalized base vector
% R     = 3x3 rotation matrix
% code  = code for signes of b and R
%           0 S(b)= UZU', R = VWU
%           1 S(b)= UZ'U', R = VWU
%           2 S(b)= UZU', R = VW'U
%           3 S(b)= UZ'U', R = VW'U
%            
% Wolfgang Förstner 8/2013
% wfoerstn@uni-bonn.de
%
% See also sugr_E_Matrix

function [b, R, code] = sugr_bR_from_E_uv(E, u, v)

% set basic matrices
W = [0 1 0; - 1, 0, 0; 0, 0, 1];
Z = [0 1 0; - 1, 0, 0; 0, 0, 0];
I = size(u, 1);

% svd of E
[U, S, V] = svd(E);
% enforce U and V to be proper rotations
U = U * det(U);
V = V * det(V);
% check sign of b and R
code = - 1;
count = 0;
for c = 0:3 % for all four cases set b and R
 
    count = count + 1;
    switch c
        case 0
            S = U * Z'*U'; R = V * W'*U';
        case 1
            S = U * Z * U';      R=V*W' * U';
        case 2
            S = U * Z'*U'; R = V * W * U';
        case 3
            S = U * Z * U';      R=V*W*U';
    end
    b = [S(3, 2); S(1, 3); S(2, 1)];
    b = b / norm(b);
    % check whether all 3D points are in direction of u and v
    sign_s = zeros(I, 1);
    sign_r = zeros(I, 1);
    for i = 1:I
        ui = u(i, :)';
        vi = v(i, :)';
        wi = R'*vi;
        m = cross(cross(b, ui), b);
        sign_s(i) = sign(det([b, m, cross(ui, wi)]));
        sign_r(i) = sign_s(i) * sign(m'*wi);
    end
    % check: the majority of points need to be in direction of u and v
    %     signs = [sign_s,sign_r];
    %     correct_sign = [mean(sign_s),mean(sign_r)];
    if mean(sign_s) > 0 && mean(sign_r) > 0
        code = c;
        return
    end
end


        
  
