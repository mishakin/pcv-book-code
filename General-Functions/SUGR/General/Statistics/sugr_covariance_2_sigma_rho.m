%% convert Covariance matrix into Standard deviation correlation matrix
%
% S_rho = sugr_covariance_2_sigma_rho(C)
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 

function S_rho = sugr_covariance_2_sigma_rho(C)

n = size(C, 1);
disp(['sigma_average: ', num2str(sqrt(trace(C) / n))]);
sigma_average = sqrt(trace(C) / n);
S_rho = zeros(n,n);
for i = 1:n
    S_rho(i, i) = sqrt(C(i, i)) / sigma_average;
    for j = 1:n
        if i ~= j
            S_rho(i, j) = C(i, j) / sqrt(C(i, i) * C(j, j));
        end
    end
end


