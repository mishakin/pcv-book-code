%% ML estimation of E-Matrix from point pairs
%
% following Sect. 13.3.5, p. 585 ff
%
% [E,sigma_0,R] = sugr_estimation_ml_E_Matrix_from_point_pairs(l,xa,T,maxiter)
%
% * l    = struct of point paris
% * xa   = 3 x 4 matrix [b,R] 
% * T    = threshold for iteration
% * maxiter = maximal iteration
%
% * E = struct estimated essential matrix, with sigma_0=1
% * sigma0 = estimated sigma0
% * R = redundancy
%
%
% Wolfgang Förstner 11/2017
% wfoerstn@uni-bonn.de
%
% See also sugr_E_Matrix, sugr_estimation_algebraic_E_Matrix_from_point_pairs

function [E, sigma_0, R] = sugr_estimation_ml_E_Matrix_from_point_pairs(l, xa, T, maxiter)

global print_option_estimation
global min_redundancy

%% Initialization
U = 5; % number of unknown parameters
Nlr = 4; % number of reduced parameters per observational group
Gc = 1; % number of constraints per observational group

lh = l.h; % Nc x Nl matrix
lCrr = l.Crr; % Nc x Nlr x Nlr matric with Crr
Nc = size(lh, 1); % number of pairs

G = Nc * Gc; % number of constraints
R = G - U; % redundancy
if R < 0
    return
end;

estl = lh; % initialize estimated observations
w_f = ones(Nc, 1); % initial weights for robust estimate

if isstruct(xa) % initiate estx, estimated unknowns
    estx = xa.bR; % [ba,Ra]
else
    estx = xa;
end

s = 0; % control variable for iterations
residuals = zeros(Nc, 1); % intial residuals

%% Start iteration -------------------------------------
for nu = 1:maxiter
    if print_option_estimation > 0
        sprintf('nu = %2', nu);
    end
    Cr = zeros(Nc, Nlr, Nlr); % reduced covariance matrices
    v_r = zeros(Nc, Nlr); % reduced residuals
    A = zeros(Nc, Gc, U); % Jacobians c -> x
    B = zeros(Nc, Gc, Nlr); % Jacobians c -> l
    W = zeros(Nc, Gc, Gc); % Weights of constraints
    cg = zeros(Nc, Gc); % constraint's residuals
 
    N_matrix = zeros(U); % normal equation matrix
    h_vector = zeros(U, 1); % right hand sides
 
    %% Build design matrices
    for n = 1:Nc
        estl_n = estl(n, :)';                 % get estl
        l_n = lh(n, :)';                   % get l
        Crr_n = squeeze(lCrr(n, :, :));
        % determine cg and Jacobians (checked)
        [lr_n, Cr_n, cg_n, atr_n, btr_n] = sugr_ghm_cg_E_matrix_from_point_pairs(l_n, estl_n, estx, Crr_n);
        % Store these
        A(n, :, :) = atr_n; % Gc x U
        B(n, :, :) = btr_n; % Gc x Nlr
        Cr(n, :, :) = Cr_n; % Nlr x Nlr
        v_r(n, :) = - lr_n';                % 1 x Nlr
        cg(n, :) = cg_n';                 % 1 x Gc
     
        % weight of contraint
        bCovb_n = btr_n * Cr_n * btr_n';      % Gc x Gc
        W(n, :, :) = inv(bCovb_n) * w_f(n); % Gc x Gc
        aW = atr_n' * squeeze(W(n,:,:)); % U x Gc
     
        N_matrix = N_matrix + aW * atr_n; % U x U
        h_vector = h_vector + aW * cg_n; % U x 1
     
    end
 
    if print_option_estimation > 0
        N_matrix = N_matrix                                                %#ok<NOPRT,ASGSL>
        l_lambda = log(eigs(N_matrix))                                     %#ok<NOPRT,NASGU>
        h_vector                                                           %#ok<NOPRT>
    end
 
    det(N_matrix);
    %% Solve
    Cxrxr = inv(N_matrix);
    estx_r = Cxrxr * h_vector;                                             %#ok<*MINV>
 
    if print_option_estimation > 0
        disp(['Result of estimation in iteration: ', num2str(nu)]);
        estimated_x = estx_r'                                               %#ok<NOPRT,NASGU>
    end
 
    max(abs(estx_r) ./ sqrt(diag(Cxrxr)));
    if max(abs(estx_r) ./ sqrt(diag(Cxrxr))) < T
        s = 2;
    end
 
    %% Determine Updates
    Omega = 0;
    % Omega_c = 0;    % check c'Wc = Omega stimmt
    check = zeros(Gc, 1);
    for n = 1:Nc
        % covariance matrix of observations (normalized)
        Clrlr = squeeze(Cr(n, :, :));
     
        % corrections of reduced observations
        delta_l_r = Clrlr * squeeze(B(n, :, :)) * squeeze(W(n, :, :)) * ...
            (cg(n, :)'-squeeze(A(n,:,:))' * estx_r) - v_r(n, :)';
        ver_r = v_r(n, :)' + delta_l_r;
     
        %sum of squared residuals
        if w_f(n) > 0
            vvp_r = ver_r' * inv(Clrlr) * ver_r;                            
            Omega = Omega + vvp_r;
            residuals(n) = vvp_r;
         
            % Omega_c = Omega_c + cg(n,:)' * squeeze(W(n,:,:)) * cg(n,:);
         
            %             % eliminate obsservation by setting w_f=0
            %             if vvp_r > 10
            %                 w_f(n)=0;
            %             end
         
        end
        % weigths = w_f'
        % [n, norm(ver_r ./ sqrt(diag(Clrlr)))];
        
        % updated estimates of both points
        estl(n, 1:3) = sugr_ghm_update_vector(estl(n, 1:3)',delta_l_r(1:2))';
        estl(n, 4:6) = sugr_ghm_update_vector(estl(n, 4:6)',delta_l_r(3:4))';
     
    end
    sigma_0 = 1;
    if R > min_redundancy
        sigma_0 = sqrt(Omega / R);
    end
    if print_option_estimation > 0
        R                                                                  %#ok<NOPRT>
        sigma_0                                                            %#ok<NOPRT>
        %sigma_c = sqrt(Omega_c/R)          % check stimmt
    end
 
    % update estimate of x
    estx = sugr_ghm_update_E_Matrix(estx, estx_r);
 
    % perform check
    for n = 1:Nc
        E = calc_S(estx(:, 1)) * estx(:, 2:4)';
        check(n) = estl(n, 1:3) * E * estl(n, 4:6)';
    end
 
    %% Stop iteration
    if s == 2
        break
    end
 
end
%% Evaluation of result ------------------------------

% determin sigma_0

if R > min_redundancy
    sigma_0 = sqrt(Omega / R);
else
    sigma_0 = 1;
end

% % choose factor
% f = 1;
% if R > min_redundancy
%     f = sigma_0;
% end

% set output
E = sugr_E_Matrix(estx(:, 1), estx(:, 2:4), Cxrxr);





