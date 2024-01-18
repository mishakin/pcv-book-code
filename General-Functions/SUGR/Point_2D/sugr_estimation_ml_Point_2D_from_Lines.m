%% ML estimate of intersection Point_2D from N Lines_2D
%
% [x,sigma_0,R] = sugr_ml_estimate_2D_Point_from_lines(l,xa,T,maxiter)
%
% * l    = struct of lines: lines, uncertain, spherically normalized homogeneous
% * xa   = struct/2-vector, approximate value
% * T    = threshold for iteration
% * maxiter = maximal iteration
%
% * x = struct estimated 2D point
% * sigma0 = estimated sigma0
% * R = redundancy
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% See also sugr_Point_2D, sugr_Line_2D, sugr_estimation_ml_Point_2D_from_Lines_Hessian
% sugr_estimation_algebraic_Point_2D_from_Lines, sugr_estimation_geometric_Point_2D_from_Lines

function [x,sigma_0,R] = ...
    sugr_estimation_ml_Point_2D_from_Lines(l,xa,T,maxiter)

global print_option_estimation
global min_redundancy

%% Initialization

lh    = l.h;
lCrr  = l.Crr;
[N,~] = size(lh);                 % number of elements
R = N-2;            % redundancy
if R < 0
    return
end;

estl = lh;                       % initialize estimated observations
w_g  = ones(N,1);                % initial weights for robust estimate

if isstruct(xa)                  % initiate estx, estimated unknowns
    estx = xa.h;
else
    estx = xa;
end

s=0;                             % control variable for iterations
residuals=zeros(N,1);            % intial residuals

%% Start iteration -------------------------------------
for nu=1:maxiter
    if print_option_estimation > 0
        sprintf('nu = %2',nu);
    end
    Cr       = zeros(N,2,2); % reduced covariance matrices
    v_r      = zeros(N,2); % reduced residuals
    A        = zeros(N,2); % Jacobians c -> x
    B        = zeros(N,2); % Jacobians c -> l
    W        = zeros(N,1); % Weights of constraints
    cg       = zeros(N,1); % constraint's residuals
    
    N_matrix = zeros(2);   % normal equation matrix
    h_vector = zeros(2,1); % right hand sides
    
    %% Build design matrices
    for n=1:N
        estl_n = estl(n,:)';                  % get estl
        l_n    =   lh(n,:)';                  % get l
        Crr_n  = squeeze(lCrr(n,:,:));
        %determine cg and Jacobians (checked)
        [lr_n,Cr_n,cg_n,atr_n,btr_n] = ...
            sugr_ghm_cg_Point_2D_from_Lines(l_n,estl_n,estx,Crr_n);
        % Store these
        A(n,:)    = atr_n(:);
        B(n,:)    = btr_n(:);
        Cr(n,:,:) = Cr_n;
        v_r(n,:)  = -lr_n';
        cg(n)     = cg_n;
        
        % weight W_g_gamma of contraint
        bCovb_n = btr_n * Cr_n * btr_n';
        W(n)    = 1 / bCovb_n * w_g(n);
        aW      = atr_n' * W(n);
        
        N_matrix = N_matrix + aW * atr_n;
        h_vector = h_vector + aW * cg_n;
        
    end
    
    det(N_matrix);
    %% Solve
    Cxrxr    = inv(N_matrix);
    estx_r   = Cxrxr * h_vector;                                           %#ok<*MINV>
    
    if print_option_estimation > 1
        disp(['Result of estimation in iteration: ',num2str(nu)]);
        %         estx_r;
    end
    
    max(abs(estx_r)./sqrt(diag(Cxrxr)));
    if max(abs(estx_r)./sqrt(diag(Cxrxr))) < T
        s=2;
    end
    
    %% Determine Updates
    Omega = 0;
    check=zeros(2,1);
    for n=1:N
        % covariance matrix of observations (normalized)
        Clrlr = squeeze(Cr(n,:,:));
        
        % corrections of reduced observations
        delta_l_r   = Clrlr * B(n,:)' * W(n) * (cg(n)-A(n,:)*estx_r)- v_r(n,:)';
        ver_r       = v_r(n,:)' + delta_l_r;
        
        % sum of squared residuals
        if w_g(n) > 0
            vvp_r = ver_r' * inv(Clrlr) * ver_r;
            Omega = Omega + vvp_r;
            residuals(n)=vvp_r;
            check=check+A(n,:)'*W(n)*B(n,:)*ver_r;
            
            % eliminate observation by setting w_g=0
            if vvp_r > 10
                w_g(n)=0;
            end
        end
        % updated estimates of observations
        estl(n,:) = sugr_ghm_update_vector(estl(n,:)',delta_l_r)';
        
    end
    if print_option_estimation > 0
        sigma_0 = sqrt(Omega/R)                                             %#ok<NOPRT,NASGU>
    end
    %     checkt = check';
    
    % update estimate of x
    %     estx0=estx;
    %     estx = sugr_ghm_update_vector(estx,estx_r);
    %     check_Atpv = (estl*estx.*w_g)';
    
    %% Stop iteration
    if s == 2
        break
    end
    
end
%% Evaluation of result ------------------------------
% determine Cxx
Cxx = null(estx') * Cxrxr * null(estx')';

% determin sigma_0
% sigma_0=1;

if R > 0
    sigma_0 = sqrt(Omega/R);
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
estCxx = f^2 * Cxx;

% set output
x = sugr_Point_2D(estx,estCxx);



