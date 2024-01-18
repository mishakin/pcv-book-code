%% ML estimate of Projection from corresponding 3D-2D point pairs
%
% [P,sigma_0,R] = sugr_estimate_ml_Projection3D_2D_from_point_pairs(X,y,xa,T,maxiter)
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
% however, 
% - internally performing conditioning/unconditioning
% - allowing for uncertain scene coordinates
%
% Wolfgang Förstner 6/2017
% wfoerstn@uni-bonn.de
%
% See also sugr_estimation_algebraic_Projection_3D_2D_from_point_pairs
% sugr_ghm_cg_Projection_3D_2D_from_point_pairs


function [P,sigma_0,R] = sugr_estimation_ml_Projection_3D_2D_from_point_pairs(X,y,xa,T,maxiter)

global print_option_estimation
global plot_option

%% Initialization
U   = 11;           % number of unknown parameters
% Nl  = 6;            % number of parameters per observational group (4+2)
Nlr = 5;            % number of reduced parameters per observational group
Gc  = 2;            % number of constraints per observational group
Nc  = size(X.e,1);  % number of pairs

%% condition points and approximate projection matrix

% condition
[Xc,MX] = sugr_condition_Points(X);
[ye,My] = sugr_condition_Points(y);
% [Xc.h,10^6*reshape(Xc.Crr,Nc,9)];
% [ye.e,10^6*reshape(ye.Cee,Nc,4)];


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
lh = zeros(Nc,6);
lCrr = zeros(Nc,5,5);
for n = 1:Nc
    lh(n,:)     = [Xc.h(n,:),ye.e(n,:)];
    lCrr(n,:,:) = [squeeze(Xc.Crr(n,:,:)) zeros(3,2); ...
        zeros(2,3) squeeze(ye.Cee(n,:,:))];
end

G      = Nc * Gc;       % number of constraints
R      = G - U;         % redundancy
if R < 0
    disp('redundancy negative')
    return
end;

nu=0;                            % Initiate iterations
estl = lh;                       % initialize estimated observations
w_f  = ones(Nc,1);               % initial weights for robust estimate

if isstruct(xa)                  % initiate estx, estimated unknowns
    estx = xa.P(:);
else
    estx = xa(:);
end


s=0;                             % control variable for iterations
residuals=zeros(Nc,Gc);          % intial residuals


%% Start iteration -------------------------------------
for nu=1:maxiter
    if print_option_estimation > 0
        sprintf('nu = %2',nu);
    end
    Cr       = zeros(Nc,Nlr,Nlr);    % reduced covariance matrices
    v_r      = zeros(Nc,Nlr);        % reduced residuals
    A        = zeros(Nc,Gc,U);       % Jacobians c -> x
    B        = zeros(Nc,Gc,Nlr);     % Jacobians c -> l
    W        = zeros(Nc,Gc,Gc);      % Weights of constraints
    cg       = zeros(Nc,Gc);         % constraint's residuals
    
    N_matrix = zeros(U);   % normal equation matrix
    h_vector = zeros(U,1); % right hand sides
    
    %% Build design matrices
    for n=1:Nc
        estl_n = estl(n,:)';                 % get estl
        l_n    = lh(n,:)';                   % get l
        Crr_n  = squeeze(lCrr(n,:,:));
        % determine cg and Jacobians (checked)
        [lr_n,Cr_n,cg_n,atr_n,btr_n] = ...
            sugr_ghm_cg_Projection_3D_2D_from_point_pairs(l_n,estl_n,estx,Crr_n);
        %cg_n'
        
        % Store these
        A(n,:,:)    = atr_n;                 % Gc x U
        B(n,:,:)    = btr_n;                 % Gc x Nlr
        Cr(n,:,:)   = Cr_n;                  % Nlr x Nlr
        v_r(n,:)    = -lr_n';                % 1 x Nlr
        cg(n,:)     = cg_n';                 % 1 x Gc
        
        % weight of contraint
        bCovb_n  = btr_n * Cr_n * btr_n';      % Gc x Gc
        W(n,:,:) = inv(bCovb_n) * w_f(n);      %#ok<MINV> % Gc x Gc                    
        aW       = atr_n' * squeeze(W(n,:,:)); % U x Gc
        
        N_matrix = N_matrix + aW * atr_n;    % U x U
        h_vector = h_vector + aW * cg_n;     % U x 1
        
    end
    
    det(N_matrix);
    log_ev=log(abs(eig(N_matrix)));
    %% Solve
    Cxrxr    = inv(N_matrix);
    estx_r   = Cxrxr * h_vector;                                           %#ok<MINV>
    
    if print_option_estimation > 0
        disp(['Result of estimation in iteration: ',num2str(nu)]);
        estimated_corr=estx_r'                                             %#ok<NOPRT,NASGU>
    end
    
%     iter_converg=[nu,max(abs(estx_r)./sqrt(diag(Cxrxr)))];
    
    if print_option_estimation > 0
        iter_maxl_minl_cond_dx=...
            [nu,max(log_ev),min(log_ev),...
            max(log_ev)-min(log_ev),max(abs(estx_r)./sqrt(diag(Cxrxr)))]   %#ok<NOPRT,NASGU>
    end
%     conv(nu)=max(abs(estx_r)./sqrt(diag(Cxrxr)));
    
    if max(abs(estx_r)./sqrt(diag(Cxrxr))) < T
        s=2;
    end
    
    %% Determine Updates
    Omega = 0;
    check = zeros(Gc,1);
    for n=1:Nc
        % covariance matrix of observations (normalized)
        Clrlr = squeeze(Cr(n,:,:));
        
        % corrections of reduced observations
        delta_l_r   = Clrlr * squeeze(B(n,:,:))' * squeeze(W(n,:,:)) ...
            * (cg(n,:)'-squeeze(A(n,:,:))*estx_r)- v_r(n,:)';
        ver_r       = v_r(n,:)' + delta_l_r;
        
        % sum of squared residuals
        if w_f(n) > 0
            vvp_r = ver_r' * inv(Clrlr) * ver_r;                           %#ok<MINV>
            Omega = Omega + vvp_r;
            residuals(n)=vvp_r;
            
            % eliminate observation by setting w_f=0 (robust version)
            %             if vvp_r > 10
            %                 w_f(n)=0;
            %             end
        end
        % updated estimates of both points
        estl(n,1:4) = sugr_ghm_update_vector(estl(n,1:4)',delta_l_r(1:3))';
        estl(n,5:6) = estl(n,5:6) + delta_l_r(4:5)';
        
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
    % check
    for n=1:Nc
        Jl = [null(estl(n,1:4))  zeros(4,2)       ; ...
            zeros(2,3)         eye(2)];
        check = check ...
            + squeeze(A(n,:,:))*null(estx')'*estx(:) ...
            + squeeze(B(n,:,:))*Jl'*estl(n,:)';
    end
%     check_cg = check';
    
    
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

