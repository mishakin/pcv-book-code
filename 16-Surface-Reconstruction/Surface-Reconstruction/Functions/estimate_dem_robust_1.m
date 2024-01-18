%% robustly estimate correction to dem, first iteration
%
% Optimization function using curvature for regularization
%   energy = sum_m rho(l_m-a_m)/sigma_n) +
%             1/sigma_d sum_ij d_ii^2+2d_ij^2+djj^2
% l_m = observed heights, interpolated bilinearly
% a_m unknown grid heights, d_ij second derivatives
%
% [A,weights,dx,v,W,Nmatrix]=...
%   estimate_dem_robust...
%    (N,U,Np,Nr,Mc,poi,xa,weights,sigma_k,iter,iter_switch,type_robust,L1,...
%    out_in,print_type,plot_type)
%
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
% type_robust = 0,1,2 (00,01,10) for points and dem
% L1          = 0/1 choice of weight function (0 = exponential, 1 = L1)
% out_in      = boolean Np-vector indicating inliers
% print_type  = 0,1,2
% plot_type   = 0,1,2
%
% A         = NxU sparse design matrix
% weights   = new weights
% dx        = delta dem
% v         = normalized residuals
% W         = weigths of obs. combining uncertainty and robust factor
% Nmatrix   = normal equation matrix
%
% wf 7/204, 4/2018
%

function [A,b,permx,weights,dx,v,W,Nmatrix]=...
    estimate_dem_robust_1...
    (N,U,Np,Nr,Mc,poi,weights,sigma_k,type_robust,...
    out_in,print_type,plot_type)

%% factor for update weights indicating outliers
k= 3;

%% start estimation
% initiate size of A
if print_type > 1
    start_time_sol = cputime;
end

A = spalloc(N,U,6*N);       % design matrix, setup once
b = zeros(N,1);             % linearized observations, zero for prior
%W = spalloc(N,N,N);
W = zeros(N,1);

if print_type > 0
    tic
end
for p = 1:Np
    x = poi(p,1);
    y = poi(p,2);
    w = 1/poi(p,4)*sqrt(weights(p))*out_in(p);
    i = floor(x);
    j = floor(y);
    aa = x-i;
    bb = y-j;
    maa = 1-aa;
    mbb=  1-bb;
    g1=maa*mbb;
    g2=aa*mbb;
    g3=maa*bb;
    g4=aa*bb;
    
    A(p,i+j*Nr+1)       = g1;
    A(p,i+j*Nr+2)       = g2;
    A(p,i+(j+1)*Nr+1)   = g3;
    A(p,i+(j+1)*Nr+2)   = g4;
    b(p)                =   poi(p,3);
    %W(p,p)=w;
    W(p)=w;
end
if print_type > 0
    display(['Time for A points  = ',num2str(toc)])
end


%% priors
if print_type > 0
    tic
end

% switch typec
%     case 0
Icc = Np;                       % first index for column curvatures
Irr = Np + (Nr-2)*Mc;                   % first index for row curvature
Irc = Np + (Nr-2)*Mc + Nr*(Mc-2);       % first index for torsion
% coefficients for column curvature
% cc  1
%    -2
%     1
Jcc = [-1,0,+1];                        % index vector
Vcc = [ 1 -2 1];                        % coefficients
wcc = 1/sigma_k;                        % /norm(Vcc)
% coefficients for row curvature
% rr  1 -2  1
Jrr = [-Nr,0,+Nr];                      % index vector
Vrr = [1 -2  1];                        % coefficients
wrr = 1/sigma_k;                        % /norm(Vrr)
% coefficients for torsion
% rc  1  -1
%    -1   1
Jrc = [-Nr-1,-Nr,-1,0];                 % indices
Vrc = [1 -1 -1 1];                      % coefficients
wrc = 1/sigma_k/sqrt(2);                % quadratic variation
% built up of priors in A
% column curvature
for j=1:Mc
    for i=2:Nr-1
        IU  = i+(j-1)*Nr;
        Icc = Icc+1;
        A(Icc,IU+Jcc) = Vcc;
        %W(Icc,Icc)    = wcc * sqrt(weights(Icc));
        W(Icc)    = wcc * sqrt(weights(Icc));
    end
end

% row curvature
for j=2:Mc-1
    for i=1:Nr
        IU  = i+(j-1)*Nr;
        Irr = Irr+1;
        A(Irr,IU+Jrr) = Vrr;
        %W(Irr,Irr)    = wrr * sqrt(weights(Irr));
        W(Irr)    = wrr * sqrt(weights(Irr));
    end
end

% torsion
for j=2:Mc
    for i=2:Nr
        IU  = i+(j-1)*Nr;
        Irc = Irc+1;
        A(Irc,IU+Jrc) = Vrc;
        %W(Irc,Irc)    = wrc * sqrt(weights(Irc));        
        W(Irc)    = wrc * sqrt(weights(Irc));
    end
end
if print_type>0
    display(['Time for A priors  = ',num2str(toc)])
end


if print_type > 2
    if U < 1000
        Ab = 2*[full(A)]
        null_A_prior = U-rank(full(A(U+1:end,:)));
        rankA = rank(full(A(U+1:end,:)));
    end  
end
%% solve

if print_type > 0
    tic
end


Aw =A;
for u=1:U
    Aw(:,u)=W.*A(:,u);
end
Nmatrix = Aw'*Aw;

if plot_type > 1
    size_N_matrix = size(Nmatrix);
    non_zeros_N_matrix =nnz(Nmatrix);
    figure
    subplot(1,3,1)
    spy(Nmatrix)
    title('N')
end

permx    = symamd(Nmatrix);
Apermw   = Aw(:,permx);

if plot_type > 1
    subplot(1,3,2)
    spy(Apermw'*Apermw)
    title('N, sorted')
end

%[C,R]    = qr(W*Aperm,W*b,0);

[C,R]    = qr(Apermw,W.*b,0);
dxperm   = R\C;
Pm       = speye(U);
Pm       = Pm(:,permx);
dx       = Pm*dxperm;

if plot_type > 1
    subplot(1,3,3)
    spy(R)
    fill_in = nnz(Rc)-(nnz(Nmatrix)/2+U/2)
    percentage_nnz   = nnz(R)/(U*(U+1)/2)
    relative_fill_in = fill_in/(nnz(Nmatrix)/2+U/2)
    title(strcat('reduced N sorted, fill in =',num2str(fill_in)))
end

if plot_type > 1
    figure
    subplot(1,2,1)
    spy(A)
    title('A')
    subplot(1,2,2)
    spy(Aperm)
    title('A sorted')
end

dx_sorted = dx;

if U < 1600 && plot_type > 0
    Aqt = R\Apermw';
    figure
    spy(Aqt')
    title('R\A sorted')
end

if print_type > 1
    time_for_solution=cputime-start_time_sol
end

% residuals (non-normalized)
v = A*dx-b;

%eso = v*(W.^2)*v)/sqrt(N-U);
eso = norm(W.*v)/sqrt(N-U);


if type_robust > 0
    
    
    % normalized residuals
    %v= W*v;
    v= W.*v;
    
    
    
    % sigma_0 for points
    eso_p = median(abs(v(1:Np)))*1.48;
    
    % sigma_0 for curvatures
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


% dem
if type_robust == 1 
    weights(Np+1:end) = min(1,1./(abs(v(Np+1:end))/kk+eps));
end

% points
if type_robust == 2 
    weights(1:Np) = min(1,1./(abs(v(1:Np))/kp+eps)).*out_in(1:Np);
end


if print_type > 1
    full(A)
end

if print_type > 0
    display(['Time for solution  = ',num2str(toc)])
end
return