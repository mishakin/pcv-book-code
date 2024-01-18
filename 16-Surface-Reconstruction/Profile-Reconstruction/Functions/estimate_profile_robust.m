%% estimate_profile_robust
%
% xest = estimate_profile_robust(N,ts,ys,sigma_e,sigma_n,
%                    Niter,type_out,type_r,print_type,plot_type)
%
% N        = scalar, length of profile
% ts,ys    = Mx1 vectors, observed profile
% sigma_e  = std of curvature
% sigma_n  = std of observations
% Niter    = number of iterations
%          = 1 -> no robust estimation
% type_out = type of outlier (0=symm, 1=asymm)
% type_r   (1) > 0 = L1, 1 = Kraus 
%          (2) = g for Kraus
%          (3) = w for Kraus
%          (4) = 1 only dem
%                2 only points
% plot_type    = 0,1,2 (no, low, high)
% print_type   = 0,1,2 (no, low, high)
%
% Wolfgang Förstner 2014-10-05
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

function [xest,A,v,weights,Cov] = estimate_profile_robust...
    (N,ts,ys,sigma_e,sigma_n,Niter,type_out,type_r,...
    print_type,plot_type)

% robust - explore sensistivity of result wrt k! Choose k in [1:3]
k = 2; 

type_robust = type_r(1);
g_factor    = type_r(2);
w_factor    = type_r(3);
rob_p_d     = type_r(4);

% Number of observations
M = length(ts);
A = zeros(M+N-2,N);

% initialize
xa = zeros(N,1);
weights = ones(M+N-2,1);
Cov = 0;
sigma_e0 = sigma_e;
gain = 1;

%% robust estimation
for iter = 1:Niter+1
    if plot_type> 0
        figure
    end
    
    sigma_e = sigma_e0*gain^(Niter-iter);
    
    if print_type > 0
        iter__sigmae_sigma_n_ratio = [iter,sigma_e,sigma_n,sigma_e/sigma_n] %#ok<NOPRT,NASGU>
        weight_before_iteration = weights(1:20)' %#ok<NOPRT,NASGU>
    end
    
    b=zeros(M+N-2,1);
    % set A and b for observed points
    for m = 1:M
        w = 1/sigma_n*sqrt(weights(m));
        A(m,ts(m)) = w;
        dl      = ys(m)-xa(ts(m));
        b(m)    = dl*w;
    end
    
    kk = M;
    % set A and b for regularizing observations
    for n = 2:N-1
        kk = kk+1;
        w = 1/sigma_e*sqrt(weights(kk));
        A(kk,n-1:n+1)= [1 -2 1]*w;
        dl           = -(xa(n+1)-2*xa(n)+xa(n-1));
        b(kk)        = dl*w;
    end
    
    % normal equation system
    Nmatrix = A'*A;
    nvector = A'*b;
    if N <= 1000
        Cov = inv(Nmatrix);
    end
    
    % solution of correction
    dx = Nmatrix\nvector;
    % corrected estimate
    xa = xa+dx;
    
%     % ----------------
%     % apply conjugate gradients
%     reorth=0;
%     k=2*N;
%     start_time=cputime
%     [dx_cg,rho,eta,F] = cgls_mod(A,b,k,reorth);
%     
%     check_conjugate_gradients = norm(dx-dx_cg)
%         time_for_solution=cputime-start_time
%     % -----------------------

    if iter == Niter+1 && print_type > 0
        xa_dx = [xa';dx'] %#ok<NOPRT,NASGU>
    end
    v = A*dx-b;                             % estimated standardized residuals
    
    %eso        = norm(v)/sqrt(M-2);              % estimated sigma0
    eso_robust = median(abs(v()))*1.48;
    
    vori=v;
    vori(1:M)      = v(1:M)*sigma_n./sqrt(weights(1:M));
    vori(M+1:end)  = v(M+1:end)*sigma_e./sqrt(weights(M+1:end));
    vs    = vori;
    
    %eso_n = median(abs(v(1:M)))*1.48;      % estimated RMS of point residuals
    eso_n = median(abs(vs(1:M)))*1.48;      % estimated RMS of point residuals   
    
    %eso_e = median(abs(v(M+1:end)))*1.48;  % estimated RMS curvature residuals
    eso_e = median(abs(vs(M+1:end)))*1.48;  % estimated RMS curvature residuals
    
    if print_type > 0
        estimated_s0_s0p_s0k=[eso,eso_p,eso_k] %#ok<NOPRT,NASGU>
    end

    % critical values
    kn = k*eso_n;
    ke = k*eso_e;

%% adapt weights
    if iter < Niter
        if type_out == 0                            % symmetric
            if print_type > 0
                disp('rob_p_d')
            end
            if iter < 4                             % first iterations
                if rob_p_d == 1 || rob_p_d == 3     % dem
                    weights(M+1:M+N-2) = ...
                        min(1,1./(abs(vs(M+1:M+N-2))/ke+0.0001));
%                     weights(M+1:M+N-2) = ...
%                         min(1,1./(abs(v(M+1:M+N-2))/ke+0.0001));
                   
                    if print_type > 1
                        iter_ke_max_v_min_w = ...
                            [iter,ke,max(abs(vs(M+1:M+N-2))),...
                            min(weights(M+1:M+N-2))]; %#ok<NASGU>
                        disp('iter_ke_max_v_min_w')
                        disp('modify dem weights')
                        residuals_new_weights_of_dem=...
                            [vs(M+1:M+N-2)';weights(M+1:M+N-2)']; %#ok<NASGU>
                        disp('residuals_new_weights_of_dem')
                    end
                    if plot_type > 0
                        subplot(3,1,1)
                        plot(2:N-1,xa(2:N-1),'-b')
                        title(num2str(sigma_e/sigma_n));
                        subplot(3,1,2);
                        hold on
                        plot(2:N-1,vs(M+1:M+N-2),'-b')
                        plot(2:N-1,vori(M+1:M+N-2),'-r')
                        plot(2:N-1,ones(N-2,1)*ke,'-b')
                        plot(2:N-1,-ones(N-2,1)*ke,'-b')
                        title(num2str(ke))
                        subplot(3,1,3);
                        plot(2:N-1,weights(M+1:M+N-2),'-b')
                        title(num2str(iter))
                    end
                end
                if rob_p_d == 2 || rob_p_d == 3     % points
                    weights(1:M) = min(1,1./(abs(vs(1:M))/kn+0.0001));
                    %weights(1:M) = min(1,1./(abs(v(1:M))/kn+0.0001));
                    
                    
                    if print_type > 1
                        residuals_new_weights_of_points_L1 = ...
                            [(1:20);v(1:20)';vs(1:20)';weights(1:20)'] %#ok<NOPRT,NASGU>
                    end
                end               
            else                                    % iteration 4 following
                if rob_p_d == 1 || rob_p_d == 3     % dem
                    weights(M+1:M+N-2) = exp(-abs(vs(M+1:M+N-2)).^2./ke^2/2)+0.0001;
                    %weights(M+1:M+N-2) = exp(-abs(v(M+1:M+N-2)).^2./ke^2/2)+0.0001;
                    
                    
                    if print_type > 1
                        iter_kk_max_v_min_w = ...
                            [iter,ke,max(abs(vs(M+1:M+N-2))),min(weights(M+1:M+N-2))] %#ok<NOPRT,NASGU>
                        disp('modify dem weights')
                        residuals_new_weights_of_dem = ...
                            [vs(M+1:M+N-2)';weights(M+1:M+N-2)'] %#ok<NOPRT,NASGU>
                    end
                    %
                    if plot_type > 0
                        subplot(3,1,1)
                        plot(2:N-1,xa(2:N-1),'-b')
                        title(num2str(sigma_e/sigma_n));
                        subplot(3,1,2);
                        hold on
                        plot(2:N-1,vs(M+1:M+N-2),'-b')
                        plot(2:N-1,vori(M+1:M+N-2),'-r')
                        plot(2:N-1,ones(N-2,1)*ke,'-b')
                        plot(2:N-1,-ones(N-2,1)*ke,'-b')
                        title(num2str(ke))
                        subplot(3,1,3);
                        plot(2:N-1,weights(M+1:M+N-2),'-b')
                        title(num2str(iter))
                    end
                end
                if rob_p_d == 2 || rob_p_d == 3         %points
                    weights(1:M) = exp(-abs(vs(1:M)).^2/kn^2/2)+0.0001;
                    %weights(1:M) = exp(-abs(v(1:M)).^2/kn^2/2)+0.0001;
                    
                    
                    if print_type > 1
                        residuals_new_weights_of_points_exp = ...
                            [(1:20);v(1:20)';vs(1:20)'/kn;weights(1:20)'] %#ok<NOPRT,NASGU>
                    end
                end                
            end
        else  %  aymmetric
            k = 1;
            g = g_factor*k;
            w = w_factor*k;
            if type_robust==0   % asymmetric L1
                weights(1:M) = max(min(1,1./abs(vs(1:M))/k+0.0001),((sign(vs(1:M))/k)+1)/2);
            else                % Kraus/Pfeifer
                for m=1:M
                    if vs(m) < g-w
                        weights(m) = 0;
                    elseif vs(m) > g
                        weights(m)=1;
                    else
                        weights(m)=1/(1+(vs(m)-g)^4);
                    end
                end
            end
        end
    else

        weights = 0.9999*ceil(weights-0.5)+0.0001;
        % initiate last iteration (Niter+1) with approximate values 0
        if iter == Niter
            xa = zeros(N,1);
        end
    end
end

xest = xa;



