%% Uncondition points
%
% x = sugr_uncondition_Points(xc,M) 
%
% xc  = struct, conditioned SUGR-points
%       xc.h = coordiante
%       xc.Crr = covariance matrices of reduced coordiantes
% M  = matrix of conditioning 
% 
% x  = struxt: unconditioned SURG-points
%
% see PCV (6.137)
% 
% Wolfgang Förstner 2/2013
% wfoerstn@uni-bonn.de 
 	
function x = sugr_uncondition_Points(xc, M)

% number and dimension
[N, d] = size(xc.h);

Mi = inv(M);

% uncondition point coordiantes
for n = 1:N
    xhcn = xc.h(n, :)';         % homogeneous condit. coord.
    Chhcn = null(xhcn') * xc.Crr(n,:,:) * null(xhcn')';
    xn = Mi * xhcn; %#ok<*MINV> % homogene ous un condit. coord.
    Chhn = Mi * Chhcn * Mi';     % ... CovM
    if d == 2
        xn = sugr_Point_2D(xn, Chhn);
    else
        xn = sugr_Point_3D(xn, Chhn);
    end
    x.h(n, :) = xn.h';        % spherically normalized uncond. coord.
    x.Crr(n, :, :) = xn.Crr; % ... reduced CovM
end

