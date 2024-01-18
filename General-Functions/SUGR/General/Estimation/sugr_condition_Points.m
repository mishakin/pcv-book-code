%% condition points, determine conditioning matrix
%
% [xc,M] = sugr_condition_Points(x)
%
% x  = struct, unconditioned points, non-homogeneous coordinates
%      x.e = coordiantes, Cartesian
%      x.Cee = covariance matrices 
%
% xc = struxt: conditioned SURG-points
%      x.h   = homogeneous coordiantes
%      x.Crr = covariance matrices of reduced coordiantes
%      x.e   = nonhomogeneous coordiante
%      x.Cee = covariance matrices of x.e
% M  = matrix of conditioning 
%
% see PCV (6.137), but sigma*sqrt(dim) instead of max
%
% Wolfgang Förstner 2/2013
% wfoerstn@uni-bonn.de 

function [xc,M] = sugr_condition_Points(x)

% number and dimension
[N,d] = size(x.e);
xc    = x;

% mean and scale
xm = mean(x.e); 
C  = cov(x.e);
sigma = sqrt(d) * sqrt(trace(C)/d);

% conditioning-matrix
M = [eye(d)       -xm';...
     zeros(1,d) sigma];
%M= eye(d+1); 

% condition point coordiantes
for n=1:N
    xhn   = [x.e(n,:) , 1]';      % homogeneous coord.
    Chhn  = [squeeze(x.Cee(n,:,:)) zeros(d,1); zeros(1,d+1)];
    xcn   = M * xhn;               % homogeneous condit. coord.
    Chhcn = M * Chhn * M';         % ... CovM
    if d ==2
        xen = sugr_Point_2D(xcn,Chhcn);
    else
        xen = sugr_Point_3D(xcn,Chhcn);
    end
    xc.h(n,:)     = xen.h';        % spherically normalized cond. coord.
    xc.Crr(n,:,:) = xen.Crr;       % ... reduced CovM
    xc.e(n,:)     = xcn(1:d)/xcn(d+1);
    xc.Cee(n,:,:) = x.Cee(n,:,:)/sigma^2;
end

%  cond_x  = cond(cov(x.e))
%  cond_xc = cond(cov(xc.h(:,1:d)./(xc.h(:,d+1)*ones(1,d))))
 

