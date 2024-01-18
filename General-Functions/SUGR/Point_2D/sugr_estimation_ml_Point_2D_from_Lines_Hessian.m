%% ML estimate of intersection Point_2D from N Lines_2D
% using line's Hessian representation
%
% model
% x cos(p) + y sin(p) - d = 0
% [x,sigma_0,R] = sugr_ml_estimate_2d_point_from_lines_Heesian(l,xa,T,maxiter)
%
% * l    = struct of lines: lines, uncertain, spherically normalized homogeneous
% * xa   = struct/2-vector, approximate value
% * T    = threshold for iteration
% * maxiter = maximal iteration
%
% * x = struct estimated point
% * sigma0 = estimated sigma0
% * R = redundancy
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% See also sugr_Point_2D, sugr_Line_2D,  sugr_estimation_ml_Point_2D_from_Lines
% sugr_estimation_algebraic_Point_2D_from_Lines, sugr_estimation_geometric_Point_2D_from_Lines

function [x,sigma_0,R] = ...
    sugr_estimation_ml_Point_2D_from_Lines_Hessian(l,xa,T,maxiter)

global print_option_estimation
global min_redundancy

%% Initialization

% select and transfer to Hessian
lh    = l.h;
lCrr  = l.Crr;

[N,~] = size(lh);                 % number of elements
R = N-2;            % redundancy
if R < 0
    return
end;

%% transfer to Hessian
Cee=zeros(N,2,2);
le = zeros(N,2);
for n=1:N
    lhn.h      = lh(n,:)';
    lhn.Crr    = reshape(lCrr(n,:,:),2,2);
    lhn.type   = 2;
    [en,Ceen]  = sugr_Line_2D_hom2Hes(lhn);
    le(n,:)    = en';
    Cee(n,1:2,1:2) = Ceen(1:2,1:2);
end

estl = le;                       % initialize estimated observations
w_g  = ones(N,1);                % initial weights for robust estimate

if isstruct(xa)                  % initiate estx, estimated unknowns
    estx = xa.h(1:2)/xa.h(3);
else
    estx = xa(1:2)/xa(3);
end

s=0;                             % control variable for iterations
residuals=zeros(N,1);            % intial residuals


%% Start iteration -------------------------------------
v       = zeros(N,2);      % Residuals
for nu=1:maxiter
    if print_option_estimation > 0
        sprintf('nu = %2',nu);
    end
    %     C        = zeros(N,2,2); % covariance matrices
    %     v_r      = zeros(N,2); % reduced residuals
    A        = zeros(N,2); % Jacobians c -> x
    B        = zeros(N,2); % Jacobians c -> l
    W        = zeros(N,1); % Weights of constraints
    cg       = zeros(N,1); % constraint's residuals
    
    N_matrix = zeros(2);   % normal equation matrix
    h_vector = zeros(2,1); % right hand sides
    
    %% Build design matrices
    for n=1:N
        estl_n = estl(n,:)';                  % get estl
        l_n    =   le(n,:)';                  % get l
        Cee_n  = squeeze(Cee(n,:,:));
        %determine cg and Jacobians (checked)
        [va,Cee_n,cg_n,at_n,bt_n] = ...
            sugr_ghm_cg_Point_2D_from_Lines_Hessian(l_n,estl_n,estx,Cee_n);
        % Store these
        A(n,:)    = at_n(:);
        B(n,:)    = bt_n(:);
        v(n,:)    = -va';
        cg(n)     = cg_n;
        
        % weight W_g_gamma of contraint
        bCovb_n = bt_n * Cee_n * bt_n';
        W(n)    = 1 / bCovb_n * w_g(n);
        aW      = at_n' * W(n);
        
        N_matrix = N_matrix + aW * at_n;
        h_vector = h_vector + aW * cg_n;
        
    end
    
    %     det(N_matrix);
    %% Solve
    Cxx    = inv(N_matrix);
    Delta_estx   = Cxx * h_vector;                                         %#ok<*MINV>
    
    if print_option_estimation > 1
        disp(['Result of estimation in iteration: ',num2str(nu)]);
        %         Delta_estx;
    end
    
    max(abs(Delta_estx)./sqrt(diag(Cxx)));
    if max(abs(Delta_estx)./sqrt(diag(Cxx))) < T
        s=2;
    end
    
    %% Determine Updates
    Omega = 0;
    check=zeros(2,1);
    for n=1:N
        % covariance matrix of observations
        Cee_n = squeeze(Cee(n,:,:));
        
        % corrections of reduced observations
        delta_l   = Cee_n * B(n,:)' * W(n) * (cg(n)-A(n,:)*Delta_estx)- v(n,:)';
        ver       = v(n,:)' + delta_l;
        
        % sum of squared residuals
        if w_g(n) > 0
            vvp = ver' * inv(Cee_n) * ver;
            Omega = Omega + vvp;
            residuals(n)=vvp;
            check=check+A(n,:)'*W(n)*B(n,:)*ver;
            
            % eliminate observation by setting w_g=0
            if vvp > 10
                w_g(n)=0;
            end
        end
        % updated estimates of observations
        estl(n,:) = estl(n,:)+delta_l';
        
    end
    if print_option_estimation > 0
        sigma_0 = sqrt(Omega/R)                                            %#ok<NOPRT,NASGU>
    end
    %     checkt = check';
    
    % update estimate of x
    %     estx0=estx;
    %     estx = estx+Delta_estx;
    %     check_Atpv = (estl*Delta_estx.*w_g)';
    
    %% Stop iteration
    if s == 2
        break
    end
    
end
%% Evaluation of result ------------------------------

% determine sigma_0
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

