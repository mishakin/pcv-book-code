%% ML estiamte of plane from points
%
% [x,sigma_0,R] = sugr_estimation_ml_symmetric_Planes_from_points(XL,XR,sigma,xa,T,maxiter)
%
% * XL, XR  = points,
% * xa      = vector, approximate value
%             xa(1:4),xa(5:8) for left and right plane
% * sigma   = std.dev. of points
% * T       = threshold for iteration
% * maxiter = maximal iteration
%
% * x         x.h estimated planes
%             x.Cxx joint covariance matrix
% * sigma0 = estimated sigma0
% * R = redundancy
%
% Wolfgang Förstner 7/2012
% wfoerstn@uni-bonn.de 

function [x, estl, sigma_0, R] = ...
    sugr_estimation_ml_symmetric_Planes_from_points(XL, XR, xa, T, maxiter)

global print_option_estimation
global min_redundancy

%% Initialization

NL = size(XL.h, 1);
NR = size(XR.h, 1);
N = NL + NR;
R = N - 4; % redundancy
if R < 0
    return
end;

lh = [XL.h; XR.h]; % initialize estimated observations
estl = lh;

lCrr = zeros(N,3,3);
for n = 1:N
    if n <= NL
        lCrr(n, :, :) = squeeze(XL.Crr(n, :, :));
    else
        lCrr(n, :, :) = squeeze(XR.Crr(n - NL, :, :));
    end
end
w_g = ones(N, 1); % initial weights for robust estimate

estx = xa;

s = 0; % control variable for iterations
residuals = zeros(N, 1); % intial residuals

%% Start iteration -------------------------------------
for nu = 1:maxiter
    if print_option_estimation > 0
        sprintf('nu = %2', nu);
    end
    Cr = zeros(N, 3, 3); % reduced covariance matrices
    v_r = zeros(N, 3); % reduced residuals
    A = zeros(N, 3); % Jacobians c -> x
    B = zeros(N, 3); % Jacobians c -> l
    W = zeros(N, 1); % Weights of constraints
    cg = zeros(N, 1); % constraint's residuals
    
    N_matrix = zeros(8); % normal equation matrix
    h_vector = zeros(8, 1); % right hand sides
    
    %% Build design matrices
    for n = 1:N
        
        estl_n = estl(n, :)';                  % get estl
        l_n = lh(n, :)';                  % get l
        Crr_n = squeeze(lCrr(n, :, :));
        %determine cg and Jacobians (checked)
        
        if n <= NL
            [lr_n, Cr_n, cg_n, atr_n, btr_n] = ...
                sugr_ghm_cg_Plane_from_points(l_n, estl_n, estx(1:4), Crr_n);
        else
            [lr_n, Cr_n, cg_n, atr_n, btr_n] = ...
                sugr_ghm_cg_Plane_from_points(l_n, estl_n, estx(5:8), Crr_n);
        end
        % Store these
        A(n, :) = atr_n(:);
        B(n, :) = btr_n(:);
        Cr(n, :, :) = Cr_n;
        v_r(n, :) = - lr_n';
        cg(n) = cg_n;
        
        % weight W_g_gamma of contraint
        bCovb_n = btr_n * Cr_n * btr_n';
        W(n) = 1 / bCovb_n * w_g(n);
        aW = atr_n' * W(n);
        if n <= NL
            N_matrix(1:3, 1:3) = N_matrix(1:3, 1:3) + aW * atr_n;
            h_vector(1:3) = h_vector(1:3) + aW * cg_n;
        else
            N_matrix(4:6, 4:6) = N_matrix(4:6, 4:6) + aW * atr_n;
            h_vector(4:6) = h_vector(4:6) + aW * cg_n;
        end
    end
    % include constraints
    factor = trace(N_matrix(1:6, 1:6));
    AL0 = estx(1:4);
    AR0 = estx(5:8);
    H = factor * ...
        [2 * [- AL0(1) * AR0(3) ^ 2, ...
        - AL0(2) * AR0(3) ^ 2, ...
        + AL0(3) * (AR0(1) ^ 2 + AR0(2) ^ 2), ...
        0] * null(AL0'),...
        2 * [+ AR0(1) * AL0(3) ^ 2, ...
        + AR0(2) * AL0(3) ^ 2, ...
        - AR0(3) * (AL0(1) ^ 2 + AL0(2) ^ 2), ...
        0] * null(AR0'),...
        [[+ AR0(2), - AR0(1), 0, 0] * null(AL0'),...
        [- AL0(2), + AL0(1), 0, 0] * null(AR0')]];
%     Hs = H / factor;
    N_matrix(7:8, 1:6) = H;
    N_matrix(1:6, 7:8) = H';
    h_vector(7:8) = factor * ...
        [- (AL0(3) ^ 2 * (AR0(1) ^ 2 + AR0(2) ^ 2) - AR0(3) ^ 2 * (AL0(1) ^ 2 + AL0(2) ^ 2)); ...
        - (AL0(1) * AR0(2) - AR0(1) * AL0(2))];
%     N_matrix / factor
%     h_vector'/factor
%     hv78 = h_vector(7:8) / factor
%     N_matrix;
%     h_vector;
    
    det(N_matrix);
    
    %% Solve
    Cxrxr = inv(N_matrix);
    estx_r = Cxrxr * h_vector;                                             %#ok<MINV>
    
    if print_option_estimation > 1
        disp(['Result of estimation in iteration: ', num2str(nu)]);
        h_vector'                                                          %#ok<NOPRT>
        estx_r'                                                            %#ok<NOPRT>
    end
    
    max(abs(estx_r(1:6)) ./ sqrt(diag(Cxrxr(1:6, 1:6))));
    if max(abs(estx_r(1:6)) ./ sqrt(diag(Cxrxr(1:6, 1:6)))) < T
        s = 2;
    end
    
    %% Determine Updates
    Omega = 0;
    check = zeros(3, 1);
    for n = 1:N
        % covariance matrix of observations (normalized)
        Clrlr = squeeze(Cr(n, :, :));
        
        % corrections of reduced observations
        if n <= NL
            delta_l_r = Clrlr * B(n, :)' * W(n) * ...
                (cg(n) - A(n, :) * estx_r(1:3)) - v_r(n, :)';
        else
            delta_l_r = Clrlr * B(n, :)' * W(n) * ...
                (cg(n) - A(n, :) * estx_r(4:6)) - v_r(n, :)';
        end
        ver_r = v_r(n, :)' + delta_l_r;
        
        % sum of squared residuals
        if w_g(n) > 0
            vvp_r = ver_r' * inv(Clrlr) * ver_r;                           %#ok<MINV>
            Omega = Omega + vvp_r;
            residuals(n) = vvp_r;
            check = check + A(n, :)'*W(n)*B(n,:)*ver_r;
            
            % eliminate observation by setting w_g=0
            if vvp_r > 10
                w_g(n) = 0;
            end
        end
        % updated estimates of observations
%         estl0 = estl;
        estl(n, :) = sugr_ghm_update_vector(estl(n, :)',delta_l_r)';
        
    end
    if print_option_estimation > 0
        sigma_0 = sqrt(Omega / R)                                          %#ok<NOPRT,NASGU>
    end
%     checkt = check';
    
    % update estimate of x
%     estx0 = estx;
    estx(1:4) = sugr_ghm_update_vector(estx(1:4), estx_r(1:3));
    estx(5:8) = sugr_ghm_update_vector(estx(5:8), estx_r(4:6));
    %check_Atpv = (estl*estx.*w_g)';
    
    %% Stop iteration
    if s == 2
        break
    end
    
end
%% Evaluation of result ------------------------------
% determine Cxx
Jrh = [null(estx(1:4)')' zeros(3, 4); zeros(3, 4) null(estx(5:8)')'];
Cxx = Jrh' * Cxrxr(1:6,1:6) * Jrh;

% determine sigma_0

if R > 0
    sigma_0 = sqrt(Omega / R);
else
    sigma_0 = 1;
end
% residuals';

% choose factor
f = 1;
if R > min_redundancy
    f = sigma_0;
end

% estimated covariance matrix of estx
estCxx = f ^ 2 * Cxx;

% estx;
% estCxx;

% set output
x.h = estx;
x.Cxx = estCxx;


