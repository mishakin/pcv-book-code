%% check estimation result, PCV Sect. 4.6.8
%
% check_estimation_result(R,x_true,S0,s0q_sample,x_sample,S,text)
%
% R          = redundancy (d.o.f. of variance factor)
% x_true     = Ux1 vector of parameters
% C0         = UxU reference covariance matrix of estimated parameters
% s0q_sample = Kx1 vector of variance factors
% x_sampl    = KxU matrix of estimated parameters
% S          = significance level of tests
% text       = string characterizing problem
%
% provides three tests
% - variance factor (if not only ones, eg for algebraic estimation)
% - covariance matrix
% - bias
%
% Wolfgang Förstner 2016-09-28
% wfoerstn@uni-bonn.de

function f = check_estimation_result(R, x_true, C0, s0q_sample, x_sample, S, text)

% number of parameters
U = length(x_true);
% sample size
K = length(find(s0q_sample(:) > 0));

disp('++++++++++++++++++++++++++++++++++++++++++++++++')
disp(strcat('Checks for ', text))
disp(['Number U of unknown parameters = ', num2str(U)]);
disp(['Redundancy R                   = ', num2str(R)]);
disp(['Number K of samples            = ', num2str(K)]);
disp('-------------------------------------------------')

%% mean and covariance matrix of sample
me = mean(x_sample);
dx = me'-x_true;
Ce = cov(x_sample);

%% check of sigma_0 (two-sided), PCV (4.354) ff
if sum(s0q_sample) ~= K
    mean_s02 = mean(s0q_sample);
    F_tol_o = finv(1 - (1 - S) / 2, K * R, 100000);
    F_tol_u = finv((1 - S) / 2, K * R, 100000);
    
    if mean_s02 > F_tol_u && mean_s02 < F_tol_o
        disp(['variance factor s_0^2      ok: mean(s_0^2) = ', ...
            num2str(mean_s02, '% 12.4f'), '     in [', ...
            num2str(F_tol_u, '% 12.4f'), ',', num2str(F_tol_o, '% 12.4f'), ']']);
    else
        disp(['variance factor  s_0^2 not ok: mean(s_0^2) = ', ...
            num2str(mean_s02, '% 12.4f'), ' not in [', ...
            num2str(F_tol_u, '% 12.4f'), ',', num2str(F_tol_o, '% 12.4f'), '] *****']);
    end
    
    % plot histogram of variance factors with expected density (F-distibution)
    screensize = plot_init;
    f = figure('Color', 'w', 'Position', [20 20 screensize / 2]);
    hold on
    % number of bins = sqrt(K)
    N_bin = ceil(sqrt(K));
    % plot histogram
    hist(s0q_sample, N_bin);
    % store bin-centres r
    [~, r] = hist(s0q_sample, N_bin);
    % determine scale
    range = abs(r(N_bin) - r(1)) * N_bin / (N_bin - 1);
    % plot density with correct scale
    plot(r, N_bin * range * fpdf(r, R, 10000), '-r', 'LineWidth', 2);
    title(['Histogram of ', num2str(K), ' variance factors with reference density (', text, ')']);
end

%% Check of covariance matrix (two-sided), PCV (4.356) ff
% Koch Parametersch Schätzung und Hypothesentests in linearen Modellen,
% 2004, (287.5)
%Ce,C0,K,S
[lambda, T] = sugr_check_CovM(Ce, C0, K, S);

if lambda > T(2) && lambda < T(1)
    disp(['covariance matrix C_xx     ok: lambda      = ', ...
        num2str(lambda, '% 12.4f'), '     in [', num2str(T(2), '% 12.4f'), ...
        ',', num2str(T(1), '% 12.4f'), ']']);
else
    disp(['covariance matrix C_xx not ok: lambda      = ', ...
        num2str(lambda, '% 12.4f'), ' not in [', num2str(T(2), '% 12.4f'), ...
        ',', num2str(T(1), '% 12.4f'), '] *****']);
end

%% Check of mean parameters (two-sided), PCV (4.359) ff

X_mean = K * dx' * inv(Ce) * dx;                                           %#ok<MINV>
X_mean_tol_h = chi2inv(1 - (1 - S) / 2, U);
X_mean_tol_l = chi2inv((1 - S) / 2, U);

if X_mean < X_mean_tol_h && X_mean > X_mean_tol_l
    disp(['mean of parameters x       ok: mean(dx)    = ', ...
        num2str(X_mean, '% 12.4f'), '     in [', ...
        num2str(X_mean_tol_l, '% 12.4f'), ',', num2str(X_mean_tol_h, '% 12.4f'), ']']);
else
    disp(['mean of parameters x   not ok: mean(dx)    = ', ...
        num2str(X_mean, '% 12.4f'), ' not in [', ...
        num2str(X_mean_tol_l, '% 12.4f'), ',', num2str(X_mean_tol_h, '% 12.4f'), ']', ' *****']);
end

