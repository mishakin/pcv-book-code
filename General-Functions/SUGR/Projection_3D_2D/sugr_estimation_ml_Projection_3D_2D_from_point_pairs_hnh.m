%% ML estimate of Projection from corresponding 3D-2D point pairs
%
% [P,sigma_0,R] = sugr_estimate_ml_Projection3D_2D_from_point_pairs_hnh(X,y,xa,T,maxiter)
%
% * X    = struct of points X.e = Nx3 matrix, X.Cee = Nx3x3 field of CovM
% * y    = struct of points y.e = Nx3 matrix, y.Cee = Nx2x2 field of CovM
% * xa   = struct/3x4-matrix, approximate value
% * T    = threshold for iteration
% * maxiter = maximal iteration
%
% * P = struct estimated homography
% * sigma0 = estimated sigma0
% * R = redundancy
%
% see Algorithm 17, p. 499,
% however internally performing conditioning/unconditioning
%
% Wolfgang Förstner 4/2018
% wfoerstn@uni-bonn.de
%
% See also sugr_estimation_algebraic_Projection_3D_2D_from_point_pairs
% sugr_ghm_cg_Projection_3D_2D_from_point_pairs


function [P,sigma_0,R] = sugr_estimation_ml_Projection_3D_2D_from_point_pairs_hnh(X,y,xa,T,maxiter)

global print_option_estimation
global plot_option

%% Initialization
U   = 11;           % number of unknown parameters
Ng  = 2;            % number of parameters per observational group
Nc  = size(X.e,1);  % number of pairs = number of constraints

%% condition points and approximate projection matrix

% condition
[Xe,MX] = sugr_condition_Points(X);
[ye,My] = sugr_condition_Points(y);


% conditioning_MX=MX
% conditioning_My=My
Pc              = condition_Homography(reshape(xa,3,4),My,MX);

% % check projection
%  xs = (Pc*Xc.h')';
%  xse = xs(:,1:2)./(xs(:,3)*ones(1,2));
%  yse = ye.e;
%  [xse-yse]

Pc = Pc/norm(Pc(:));
xa = Pc(:);

%% provide data for estimation
le = zeros(Nc,2);
lCee = zeros(Nc,2,2);
for n = 1:Nc
    le(n,:)     = ye.e(n,:);
    lCee(n,:,:) = squeeze(ye.Cee(n,:,:));
end

N      = Nc * Ng;       % number of observations
R      = N - U;         % redundancy
if R < 0
    disp('redundancy negative')
    return
end;

nu=0;                            % Initiate iterations
estl = le;                       % initialize estimated observations

if isstruct(xa)                  % initiate estx, estimated unknowns
    estx = xa.P(:);
else
    estx = xa(:);
end


s=0;                             % control variable for iterations
residuals=zeros(Nc,Ng);          % intial residuals


%% Start iteration -------------------------------------
for nu=1:maxiter
    if print_option_estimation > 0
        sprintf('nu = %2',nu);
    end
    C        = zeros(Nc,Ng,Ng);      % covariance matrices
    v        = zeros(Nc,Ng);         % residuals
    A        = zeros(Nc,Ng,U);       % Jacobians c -> x
    W        = zeros(Nc,Ng,Ng);      % Weights of constraints
    dl       = zeros(Nc,Ng);         % linearized observations
    
    N_matrix = zeros(U);   % normal equation matrix
    h_vector = zeros(U,1); % right hand sides
    
    %% Build design matrices
    % 12x11 Jacobian for xa
    Jr     = null(xa');
    
    for n=1:Nc
        l_n    = le(n,:)';                   % get l
        C_n    = squeeze(lCee(n,:,:));       % get CovM(l)
        % approximate fitted point homogeneous/non-homogeneous
        est_lh_n   = reshape(estx,3,4)*[Xe.e(n,:)';1];
        est_l_n    = est_lh_n(1:2)/est_lh_n(3);
        % determine dl
        dl_n   = l_n - est_l_n;       
        % 2x3 Jacobian for c(x)
        Jc_n     = [est_lh_n(3)*eye(2), -est_lh_n(1:2)]/est_lh_n(3)^2;
        % 2x11 design matrix for n-th observation 
        atr_n  = Jc_n*kron([Xe.e(n,:),1],eye(3))*Jr;

        % Store these
        A(n,:,:)    = atr_n;                 % Ng x U
        C(n,:,:)    = C_n;                   % Ng x Ng
        dl(n,:)     = dl_n';                 %  1 x Ng
        
        % weight of contraint
        W(n,:,:) = inv(C_n);                   % Nc x Nc                  
        aW       = atr_n' * squeeze(W(n,:,:)); %  U x Ng
        
        % normal equation system
        N_matrix = N_matrix + aW * atr_n;    % U x U
        h_vector = h_vector + aW * dl_n;     % U x 1
        
    end
    
            det(N_matrix);
            log_ev=log(abs(eig(N_matrix)));
    %% Solve for reduced parameters of the projection matrix
    Cxrxr    = inv(N_matrix);
    estx_r   = Cxrxr * h_vector;                                           %#ok<MINV>
    
            if print_option_estimation > 0
                disp(['Result of estimation in iteration: ',num2str(nu)]);
                estimated_corr=estx_r'  
                iter_maxl_minl_cond_dx=...
                    [nu,max(log_ev),min(log_ev),...
                    max(log_ev)-min(log_ev),max(abs(estx_r)./sqrt(diag(Cxrxr)))]   %#ok<NOPRT,NASGU>
            end
    % check for convergence
    if max(abs(estx_r)./sqrt(diag(Cxrxr))) < T
        s=2;
    end
    
    %% Determine Updates
    Omega = 0;
    for n=1:Nc
        % covariance matrix of observations (normalized)
        C_n = squeeze(C(n,:,:));
        % fitted observation
        est_lh_n   = reshape(estx,3,4)*[Xe.e(n,:)';1];
        est_le_n   = est_lh_n(1:2)/est_lh_n(3);
        % residuals
        ver        = est_le_n - le(n,:)';
        % sum of squared residuals
        vvp = ver' * inv(C_n) * ver;                           %#ok<MINV>
        Omega = Omega + vvp;
        residuals(n)=vvp;
    end
    
    sigma_0 = 1;
    if R > 0
        sigma_0 = sqrt(Omega/R);
    end
    if print_option_estimation > 0
        sigma_0                                                            %#ok<NOPRT>
    end
    
    % update estimate of x
%     estx0 = estx;
    estx  = sugr_ghm_update_vector(estx,estx_r);
    %estx'
    % check A' W v = 0   
    check = zeros(U,1);
    for n=1:Nc
        check = check ...
            + squeeze(A(n,:,:))'*inv( squeeze(C(n,:,:)))*dl(n,:)';
    end
    check_cg = check';
    
    
    %% Stop iteration
    if s == 2
        break
    end
    
end

if plot_option == -1
    figure(2)
    hold on
    plot(1:nu,log(conv),'bo-')
end
%% Evaluation of result ------------------------------

% log_condition_number_conditioned=...
%     max(log_ev)-min(log_ev);

% determine Cxx
Jrh = null(estx');
Cxx= Jrh * Cxrxr * Jrh';                                                   %#ok<MINV>

% determine sigma_0

if R > 0
    sigma_0 = sqrt(Omega/R);
else
    sigma_0 = 1;
end

% set output, still conditioned
Pce = reshape(estx,3,4);
Pc  = sugr_Projection_3D_2D(Pce,Cxx);

%% uncondition projection matrix
P = sugr_uncondition_Projection(Pc,My,MX);

