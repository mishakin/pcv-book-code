%% generate_observed_AR2 with noise and ouliers
%
% [x,y,out_in,select,xs,ys] = generate_observed_AR2(N, sigma_e, sigma_n, Pout, Mout, type_out, dens)
%
% N         = number of points
% sigma_e   = process noise
% sigma_n   = observation noise
% Pout      = probability for outliers
% Mout      = maximal outlier
% type_out  = type otliers (0=symm., 1=asymm.)
% dens      = percentage of observed points in (0,1]
%
% x         = true profile with slope 0
% y         = noisy profile with slope 0 (no outliers)
% out_in    = Nx1 vector, 0/1: o inlier, 1 outlier
% select    = index of observed points
% xs        = observed true profile
% ys        = observed noisy profile with outliers
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 

function [x, y, out_in, select, xs, ys] = ...
    generate_observed_AR2(N, sigma_e, sigma_n, Pout, Mout, type_out, dens)

% initiate true profiles
x = zeros(N, 1);
y = zeros(N, 1);

% generate true and noisy profile
for n = 3:N
    x(n) = 1.9998 * x(n - 1) - 0.9999 * x(n - 2) + randn(1) * sigma_e;
    y(n) = x(n) + randn(1) * sigma_n;
end

% enforce slope 0
for i = 1:N
    y(i) = y(i) - (i - 1) * (x(N) - x(1)) / (N - 1);
    x(i) = x(i) - (i - 1) * (x(N) - x(1)) / (N - 1);
end
% add ouliers and select observations
out_in = zeros(N, 1);
m = 0;
for i = 1:N
    % add outliers
    if rand(1) < Pout
        out_in(i) = 1;
        if type_out == 0 % symmetric in [-Mout,+Mout]
            y(i) = y(i) + 2 * (rand(1) - 0.5) * Mout;
        else % asymmetric in [0,Mout]
            y(i) = y(i) + rand(1) * Mout;
        end
    end
    
    % select observations
    if rand(1) < dens
        m = m + 1;
        select(m) = i;                                                      %#ok<*AGROW>
        xs(m) = x(i);
        ys(m) = y(i);
    end
end



