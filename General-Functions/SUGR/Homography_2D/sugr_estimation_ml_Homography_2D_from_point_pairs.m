%% ML estimate of 2D homography from point pairs
%
% [H,sigma_0,R] = sugr_estimate_ml_Homography_2D_from_point_pairs(l,xa,T,maxiter)
%
% * l    = struct of point paris
% * xa   = struct/3x3-matrix, approximate value
% * T    = threshold for iteration
% * maxiter = maximal iteration
%
% * H = struct estimated homography
% * sigma0 = estimated sigma0
% * R = redundancy
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%

function [H, sigma_0, R] = sugr_estimation_ml_Homography_2D_from_point_pairs(l, xa, T, maxiter)

global print_option_estimation

%% Initialization
U = 8; % number of unknown parameters
Nlr = 4; % number of reduced parameters per observational group
Gc = 2; % number of constraints per observational group

lh = l.h;
lCrr = l.Crr;
Nc = size(lh, 1); % number of pairs

G = Nc * Gc; % number of constraints
R = G - U; % redundancy
if R < 0
    return
end

estl = lh; % initialize estimated observations
w_f = ones(Nc, 1); % initial weights for robust estimate

if isstruct(xa) % initiate estx, estimated unknowns
    estx = xa.H;
else
    estx = xa;
end

s = 0; % control variable for iterations
residuals = zeros(Nc, 2); % intial residuals

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
        [lr_n, Cr_n, cg_n, atr_n, btr_n] = sugr_ghm_cg_2D_homography_from_point_pairs(l_n, estl_n, estx, Crr_n);
        % Store these
        A(n, :, :) = atr_n; % Gc x U
        B(n, :, :) = btr_n; % Gc x Nlr
        Cr(n, :, :) = Cr_n; % Nlr x Nlr
        v_r(n, :) = - lr_n';                % 1 x Nlr
        cg(n, :) = cg_n';                 % 1 x Gc
        
        % weight of contraint
        bCovb_n = btr_n * Cr_n * btr_n';      % Gc x Gc
        W(n, :, :) = inv(bCovb_n) * w_f(n); %#ok<MINV> % Gc x Gc
        aW = atr_n' * squeeze(W(n,:,:)); % U x Gc
        
        N_matrix = N_matrix + aW * atr_n; % U x U
        h_vector = h_vector + aW * cg_n; % U x 1
        
    end
    
    %     det(N_matrix);
    %% Solve
    Cxrxr = inv(N_matrix);
    estx_r = Cxrxr * h_vector; %#ok<MINV>
    
    if print_option_estimation > 0
        disp(['Result of estimation in iteration: ', num2str(nu)]);
        disp(estx_r);
    end
    
    %     max(abs(estx_r)./sqrt(diag(Cxrxr)));
    if max(abs(estx_r) ./ sqrt(diag(Cxrxr))) < T
        s = 2;
    end
    
    %% Determine Updates
    Omega = 0;
    check = zeros(Gc, 1);
    for n = 1:Nc
        % covariance matrix of observations (normalized)
        Clrlr = squeeze(Cr(n, :, :));
        
        % corrections of reduced observations
        delta_l_r = Clrlr * squeeze(B(n, :, :))' * squeeze(W(n,:,:)) * (cg(n,:)' - squeeze(A(n, :, :)) * estx_r) - v_r(n, :)';
        ver_r = v_r(n, :)' + delta_l_r;
        
        % sum of squared residuals
        if w_f(n) > 0
            vvp_r = ver_r' * inv(Clrlr) * ver_r;                           %#ok<MINV>
            Omega = Omega + vvp_r;
            residuals(n) = vvp_r;
            
            % eliminate obsservation by setting w_f=0
            if vvp_r > 10
                w_f(n) = 0;
            end
        end
        % updated estimates of both points
        estl(n, 1:3) = sugr_ghm_update_vector(estl(n, 1:3)',delta_l_r(1:2))';
        estl(n, 4:6) = sugr_ghm_update_vector(estl(n, 4:6)',delta_l_r(3:4))';
        
    end
    sigma_0 = 1;
    if R > 0
        sigma_0 = sqrt(Omega / R);
    end
    if print_option_estimation > 0
        disp(sigma_0)
    end
    if print_option_estimation > 1
        disp(check')
    end
    
    % update estimate of x
    if print_option_estimation > 1
        disp(estx)
    end
    
    estx = sugr_ghm_update_homography(estx, estx_r);
    for n = 1:Nc
        Jl = [null(estl(n, 1:3)) zeros(3, 2); zeros(3, 2) null(estl(n, 4:6))];
        [~, Jhr] = sugr_minimal_Homography_2D(estx, zeros(9));
        check = check ...
            + squeeze(A(n, :, :)) * Jhr * estx(:) ...
            + squeeze(B(n, :, :)) * Jl'*estl(n,:)';
    end
    if print_option_estimation > 1
        disp(check)
    end
    
    %% Stop iteration
    if s == 2
        break
    end
    
end

%% Evaluation of result ------------------------------
% determine Cxx

[~, Jrh] = sugr_get_CovM_homogeneous_Homography_2D(sugr_Homography_2D(estx, zeros(9)));
Cxx = Jrh * Cxrxr * Jrh';                                                  %#ok<MINV>

% determine sigma_0
if R > 0
    sigma_0 = sqrt(Omega / R);
else
    sigma_0 = 1;
end
if print_option_estimation > 1
    disp(residuals');
end

% % choose factor
% f = 1;
% if R > min_redundancy
%     f = sigma_0;
% end

% set output
H = sugr_Homography_2D(estx, Cxx);


