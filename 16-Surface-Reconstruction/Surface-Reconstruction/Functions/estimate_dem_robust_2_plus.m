%% robustly estimate correction to dem, using design matrix
%
% Optimization function using curvature for regularization
%   energy = sum_m rho(l_m-a_m)/sigma_n) +
%             1/sigma_d sum_ij d_ii^2+2d_ij^2+djj^2
% l_m = observed heights, interpolated bilinearly
% a_m unknown grid heights, d_ij second derivatives
%
% [weights,dx,v,W,Nmatrix]=...
%     estimate_dem_robust_2_plus...
%     (A,b,N,U,Np,Nr,Mc,poi,xa,weights,sigma_k,iter,iter_switch,type_robust,L1,...
%     out_in,print_type,plot_type)
%
% A           = design matrix (sparse)
% b           = right hand sides
% N           = number of observations (points+second derivatives)
% U           = number of unknowns (dem)
% Np          = number of points
% Nr          = number of rows
% Mc          = number oc columns
% poi         = Npx4 matrix of points (x,y,h,sigma), x,y, referring to grid
% weights     = Nx1 vector of weights in [0,1]
% sigma_k     = standard deviation of second derivatives
% iter        = current iteration for controling weighting
% iter_switch = iteration number where to switch type of weighting
% type_robust = 0,1,2, (00,01,10) for points and dem
% L1          = 0/1 choice of weight function (0 = exponential, 1 = L1)
% out_in      = boolean Np-vector indicating inliers
% print_type  = 0,1,2,3
% plot_type   = 0,1,2
%
% weights   = new weights
% dx        = delta dem
% v         = normalized residuals
% W         = weigths of obs. combining uncertainty and robust factor
% Nmatrix   = normal equation matrix
%
% wf 7/204, 4/2018
%

function [weights,dx,v,W,Nmatrix]=...
    estimate_dem_robust_2_plus...
    (A,b,permx,N,U,Np,Nr,Mc,poi,weights,sigma_k,type_robust,L1,...
    out_in,print_type)

Nmatrix=0;

%% factor for update weights indicating outliers
k= 3;

%% start estimation
if print_type > 1
    start_time_sol = cputime;
end

%b = zeros(N,1);
%W = spalloc(N,N,N);
W = zeros(N,1);

% if print_type>0
%     tic
% end

%% determine linearized observations dl by bilinear intoerpolation
for p = 1:Np
   %W(p,p)= 1/poi(p,4)*sqrt(weights(p))*out_in(p);
   W(p)= 1/poi(p,4)*sqrt(weights(p))*out_in(p);
end

% if print_type > 0
%     display(['Time for dl points = ',num2str(toc)])
% end

%% determin lineraized observations for priors

% if print_type > 0
%     tic
% end
% switch typec
%     case 0
Icc = Np;                            % first index for column curvatures
Irr = Np + (Nr-2)*Mc;                % first index for row curvature
Irc = Np + (Nr-2)*Mc + Nr*(Mc-2);    % first index for torsion
% coefficients for column curvature
% cc  1
%    -2
%     1
Jcc = [-1,0,+1];                     % index vector
Vcc = [ 1 -2 1];                     % coefficients
wcc = 1/sigma_k;                     %/sqrt(6);  % /norm(Vcc)
% coefficients for row curvature
% rr  1 -2  1
Jrr = [-Nr,0,+Nr];                   % index vector
Vrr = [1 -2  1];                     % coefficients
wrr = 1/sigma_k;                     % /norm(Vrr)
% coefficients for torsion
% rc  1  -1
%    -1   1 ^
Jrc = [-Nr-1,-Nr,-1,0];              % indices
Vrc = [1 -1 -1 1];                   % coefficients
wrc = 1/sigma_k/sqrt(2);             % for quadratic variation
% built up of priors in A
% column curvature
for j=1:Mc
    for i=2:Nr-1
        Icc = Icc+1;
       %W(Icc,Icc)    = wcc * sqrt(weights(Icc));
        W(Icc)    = wcc * sqrt(weights(Icc));
    end
end
% row curvature
for j=2:Mc-1
    for i=1:Nr
        Irr = Irr+1;
        %W(Irr,Irr)    = wrr * sqrt(weights(Irr));
        W(Irr)    = wrr * sqrt(weights(Irr));
    end
end
% torsion
for j=2:Mc
    for i=2:Nr
        Irc = Irc+1;
        %W(Irc,Irc)    = wrc * sqrt(weights(Irc));
        W(Irc)    = wrc * sqrt(weights(Irc));
    end
end
% if print_type > 0
%     display(['Time for dl priors = ',num2str(toc)])
% end

if print_type > 1
    time_for_building_dl=toc
end

%% solve

if print_type > 0
    tic
end

Aperm    = A(:,permx);

%[C,R]    = qr(W*A,W*b,0);

Apermw =Aperm;
for u=1:U
    Apermw(:,u)=W.*Aperm(:,u);
end

[C,R]    = qr(Apermw,W.*b,0);
dxperm   = R\C;
Pm       = speye(U);
Pm       = Pm(:,permx);
dx       = Pm*dxperm;

% residuals (non-normalized)
v = A*dx-b;

% estimated sqrt of variance factor
eso = norm(W.*v)/sqrt(N-U);

if type_robust > 0
    % robust sigma_0 for points
    eso_p = median(abs(v(1:Np)))*1.48;
    
    % robust sigma_0 for curvatures
    eso_k = median(abs(v(Np+1:end)))*1.48;
    if print_type > 1
        est_s0_s0p_s0k=[eso,eso_p,eso_k]
    end
     
    % critical values
    kp = k*eso_p;
    kk = k*eso_k;
      
    if print_type > 1
        criticalvalue_kp_sigmap_kk_sigma_k=[kp,poi(1,4),kk,sigma_k]
    end
end


%% perform reqeighting on points or dem
if L1  % for the first iterations
    % dem
    if type_robust == 1
        weights(Np+1:end) = min(1,1./(abs(v(Np+1:end))/kk+eps));
    end
    % points
    if type_robust == 2
        weights(1:Np) = min(1,1./(abs(v(1:Np))/kp+eps)).*out_in(1:Np);
    end
else
    if type_robust == 1
        weights(Np+1:end) = exp(-abs(v(Np+1:end)).^2/kk^2/2)+0.0001;
    end
    if type_robust == 2
        weights(1:Np) = (exp(-abs(v(1:Np)).^2/kp^2/2)+0.0001).*out_in(1:Np);
    end
end


if print_type > 0
    display(['Time for solution  = ',num2str(toc)])
end

return