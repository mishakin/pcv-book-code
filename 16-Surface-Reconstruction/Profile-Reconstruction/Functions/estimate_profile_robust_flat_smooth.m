%% estimate_profile_robust_flat_smooth
%
% xest = estimate_profile_robust_flat_smooth(N,fs,ts,ys,
%          sigma_e1,sigma_e2,sigma_n,Niter,type_out,type_r,print_type)
%
% N        = scalar, length of profile
% fs       = Nx1 vector 1,2
% ts,ys    = Mx1 vectors, observed profile
% sigma_e  = std of curvature
% sigma_n  = std of observations
% Niter    = number of iterations
% type_out = type of outlier (0=symm, 1=asymm)
% type_r   (1) > 0 = L1, 1 = Kraus 
%          (2) = g for Kraus
%          (3) = w for Kraus
%          (4) = 1 only dem
%                2 only points
%                3 both
% print_type   = 0,1,2 (no, low, high)
%
% Wolfgang Förstner 2014-10-05
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

function [xest,A,v,weights,Cov] = estimate_profile_robust_flat_smooth(...
    N,fs,ts,ys,sigma_e1,sigma_e2,sigma_n,Niter,type_out,type_r,print_type)

% robust
type_robust=type_r(1);
g_factor   =type_r(2);
w_factor   =type_r(3);
rob_p_d    =type_r(4);

% Number of observations
M = length(ts);
A=zeros(M+N-2,N);

% initialize
xa=zeros(N,1);
weights=ones(M+N-2,1);
Cov=0;

%% robust estimation
for iter = 1:Niter+1
    
    if print_type > 0
        disp(iter)
    end
    
    b = zeros(M+N-1,1);
    for m = 1:M
        w = 1/sigma_n*sqrt(weights(m));
        A(m,ts(m)) = w;
        dl      = ys(m)-xa(ts(m));
        b(m)    = dl*w;
    end
    
    if iter == Niter+1 && print_type > 1
        rhs_w = [b(1:M)';weights(1:M)'];
        disp(rhs_w)
    end
    
    k = M;
    for n = 1:N-1
        k = k+1;
        if n == 1 || fs(n) == 1 
            A(k,n:n+1)   = [1 -1]/sigma_e1;
        else
            A(k,n-1:n+1) = [1 -2 1]/sigma_e2;
        end
    end
    
    % normal system
    Nmatrix = A'*A;
    nvector = A'*b;
    if N <= 1000 
        Cov = inv(Nmatrix);
    end
    dx = Nmatrix\nvector;
    xa = xa+dx;
    
    
    if print_type > 0
        norm_dx=[iter,norm(dx)];
        disp(norm_dx)
    end
    
    if iter == Niter+1 && print_type > 1
        xa_dx=[xa';dx'];
        disp(xa_dx)        
    end
    v = A*dx-b;
    
    %% adapt weights
    if print_type > 0
        disp('adapt weights')
        disp([iter,rob_p_d])
        sv = median(abs(v(M+1:M+N-2)))*1.48;
        disp(sv)
    end
    
    if iter < Niter

        if type_out == 0
        
            k=3;
            if iter < 4
                if rob_p_d == 2 || rob_p_d == 3
                    weights(1:M) = min(1,1./abs(v(1:M))/k+0.0001);
                else
                    weights(M+1:M+N-2) = min(1,1./abs(v(M+1:M+N-2))/k+0.0001);
                end
            else
                if rob_p_d == 2 || rob_p_d == 3
                    weights(1:M) = exp(-abs(v(1:M)).^2/k^2);
                else
                    weights(M+1:M+N-2) = exp(-abs(v(M+1:M+N-2)).^2/(k^2*sv^2));
                end
            end
            
        else  % asymmetric L1 or Kraus/Pfeifer  
            
            k = 1;
            g = g_factor*k;
            w = w_factor*k;
            
            if type_robust == 0
                weights(1:M) = max(min(1,1./abs(v(1:M))/k+0.0001),((sign(v(1:M))/k)+1)/2);
                
            else
                for m = 1:M
                    if v(m) < g-w
                        weights(m) = 0;
                        
                    elseif v(m) > k
                        weights(m)=1;
                        
                    else
                        weights(m)=1/(1+(v(m)-g)^4);
                        
                    end
                end
            end
        end
    else
        
        weights(1:M) = weights(1:M)>0.5;
        if iter == Niter
            xa = zeros(N,1);
        end
    end
    
end

xest = xa;


