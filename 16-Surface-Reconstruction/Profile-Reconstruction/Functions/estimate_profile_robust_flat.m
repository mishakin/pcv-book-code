%% estimate_profile_robust_flat
%
% xest = estimate_profile_robust(N,ts,ys,sigma_e,sigma_n,
%                    Niter,type_out,type_r,print_type)
%
% N        = scalar, length of profile
% ts,ys    = Mx1 vectors, observed profile
% sigma_e  = std of curvature
% sigma_n  = std of observations
% Niter    = number of iterations
% type_out = type of outlier (0=symm, 1=asymm)
% type_r   (1) > 0 = L1, 1 = Kraus 
%          (2) = g for Kraus
%          (3) = w for Kraus
% print_type   = 0,1,2 (no, low, high)
%
% Wolfgang Förstner 2014-10-05
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

function [xest,A,v,weights,Cov] = estimate_profile_robust_flat...
    (N,ts,ys,sigma_e,sigma_n,Niter,type_out,type_r,print_type)

% robust
type_robust=type_r(1);
g_factor   =type_r(2);
w_factor   =type_r(3);

% Number of observations
M = length(ts);
A = zeros(M+N-2,N);

% initialize
xa = zeros(N,1);
weights = ones(M+N-2,1);
Cov = 0;

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

    k = M;
    for n = 1:N-1
        k = k+1;
        A(k,n:n+1)=[1 -1]/sigma_e;
    end
    
    % normal system
    Nmatrix = A'*A;
    nvector = A'*b;
    if N <= 1000 
        Cov = inv(Nmatrix);
    end
    dx = Nmatrix\nvector;
    xa = xa+dx;
    
    if iter == Niter+1 && print_type > 0
        xa_dx = [xa';dx'];
        disp(xa_dx)
    end
    v = A*dx-b;
    
    %% adapt weights
    if iter < Niter
        
        if type_out == 0
            
            k=1;
            if iter < 6
                weights(1:M) = min(1,1./abs(v(1:M))/k+0.0001);
            else
                weights(1:M) = exp(-abs(v(1:M))/k);
            end
            
        else  % Kraus/Pfeifer  
            
            k = 1;
            g = g_factor*k;
            w = w_factor*k;
            
            if type_robust == 0
                
                weights(1:M) = ...
                    max( min(1,1./abs(v(1:M))/k+0.0001),...
                    ((sign(v(1:M))/k)+1)/2 );
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
            xa  = zeros(N,1);
        end
    end
    
end

xest = xa;


